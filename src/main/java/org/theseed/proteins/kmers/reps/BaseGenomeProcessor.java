package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.IOFileFilter;
import org.apache.commons.io.filefilter.RegexFileFilter;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.TextStringBuilder;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.kmers.IRepGenContainer;
import org.theseed.proteins.kmers.ProteinData;
import org.theseed.proteins.kmers.ProteinDataFactory;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.stats.GenomeEval;
import org.theseed.stats.QualityCountMap;

/**
 * This is the base class for both GenomeProcessor and UpdateProcessor.  It performs the common tasks
 * in maintaining the various files that depend on genome evaluations.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseGenomeProcessor extends BaseProcessor implements IRepGenContainer {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(BaseGenomeProcessor.class);
    /** manager for key data about each genome and its seed protein */
    private ProteinDataFactory genomeList;
    /** list of repgen sets */
    private List<RepGenomeDb> repGenSets;
    /** counters for the size of each repgen group */
    private List<QualityCountMap<String>> statMaps;
    /** repGen set used for close-genome finder */
    private RepGenomeDb lastRepGen;
    /** groups for last repgen set */
    private Map<String, List<ProteinData>> repFinderSets;
    /** map of representative genomes to the directories where they belong */
    private final Map<String, List<File>> repHash;
    /** header line for four-column table containing SSUs */
    private static final String FOUR_COLUMN_HEADER = "genome_id\tgenome_name\tseed_protein\tssu_rna\tseed_dna";
    /** file name pattern for GTOs */
    private static final Pattern GTO_FILE_PATTERN = Pattern.compile("\\d+\\.\\d+\\.gto");
    /** file name pattern for repgen GTO directories */
    private static final Pattern GTO_DIR_PATTERN = Pattern.compile("GTO\\d+");
    /** file name filter for GTO files */
    private static final IOFileFilter GTO_FILE_FILTER = new RegexFileFilter(GTO_FILE_PATTERN);
    /** file name filter for GTO directories */
    private static final IOFileFilter GTO_DIR_FILTER =  new RegexFileFilter(GTO_DIR_PATTERN);

    public enum GtoScheme {
        /** download GTOs */
        YES,
        /** do not download GTOs */
        NO,
        /** do not build repgens; only download GTOs */
        ONLY;
    }

    // COMMAND-LINE OPTIONS

    /** specifies whether bad-SSU genomes will be considered good */
    @Option(name = "--minLevel", usage = "minimum rating level to keep")
    private ProteinData.Rating minRating;

    /** number of genomes per batch when retrieving sequences */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "200", usage = "number of genomes per batch")
    private int batchSize;

    /** processing mode, GTOs only, everything, no GTOs */
    @Option(name = "--gtos", usage = "processing mode relating to GTO download")
    private GtoScheme gtoMode;

    /** output directory */
    @Argument(index = 0, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    /**
     * @return the list of representative-genome databases
     */
    protected List<RepGenomeDb> getRepGenSets() {
        return repGenSets;
    }

    /**
     * @return the batch size
     */
    protected int getBatchSize() {
        return batchSize;
    }

    /**
     * @return the output directory
     */
    protected File getOutDir() {
        return outDir;
    }

    public BaseGenomeProcessor() {
        super();
        this.repHash = new TreeMap<>();
    }

    /**
     * Set the base-class option defaults.
     */
    protected void setBaseDefaults() {
        this.batchSize = 500;
        this.minRating = ProteinData.Rating.SINGLE_SSU;
        this.gtoMode = GtoScheme.YES;
    }

    /**
     * Separate the genomes into RepGen sets.  A set of genomes that do not need to be checked
     * can be specified to improve performance in updates. This is the place where the representative
     * genomes are chosen.
     * 
     * We do one pass for each repgen set, starting with the smallest. Each larger set starts with
     * the previous set's genomes pre-selected, so that the smaller set is always a subset of the
     * larger one. To make this efficient, 
     *
     * @param skipSet	IDs of genomes that do not need to be checked
     */
    protected void collateGenomes(Set<String> skipSet) {
        log.info("Sorting genomes into repgen sets.");
        Iterator<RepGenomeDb> iter = this.repGenSets.iterator();
        RepGenomeDb repGenDb = iter.next();
        this.buildRepGenSet(repGenDb, skipSet);
        while (iter.hasNext()) {
            RepGenomeDb newRepGenDb = iter.next();
            copyReps(repGenDb, newRepGenDb);
            repGenDb = newRepGenDb;
            // We can't skip anything for the larger sets. Because of the subsetting requirement,
            // we only get to reuse the data from the smallest set.
            this.buildRepGenSet(repGenDb, Collections.emptySet());
        }
    }

    /**
     * Copy all the representatives from the incoming repgen set to a new repgen database.
     * 
     * @param repGenDb      source repgen database
     * @param newRepGenDb   target repgen database
     */
    public static void copyReps(RepGenomeDb repGenDb, RepGenomeDb newRepGenDb) {
        for (RepGenome rep : repGenDb)
            newRepGenDb.addRep(rep);
        log.info("{} representatives copied from {} into {}.", repGenDb.size(), repGenDb, newRepGenDb);
    }

    /**
     * Check all the genomes not in the skip set to find out if they qualify as representatives.
     * 
     * @param repGenDb  target representative-genome database
     * @param skipSet   set of genome IDs for genomes to skip
     */
    private void buildRepGenSet(RepGenomeDb repGenDb, Set<String> skipSet) {
        long start = System.currentTimeMillis();
        int batchCount = this.batchSize;
        int sim = repGenDb.getThreshold();
        int genomeCount = 0;
        // First, make all genomes in the set represent themselves. This gives us a minor performance
        // boost.
        log.info("Initializing default representation for rep{}.", sim);
        for (RepGenome rep : repGenDb) {
            String repId = rep.getGenomeId();
            ProteinData genomeData = this.genomeList.getGenome(repId);
            genomeData.setRepresentation(repGenDb, repId, RepGenome.INFINITY, 0.0);
        }
        log.info("Processing rep{} set with {} starting representatives.", sim, repGenDb.size());
        for (ProteinData genomeData : this.genomeList) {
            genomeCount++;
            if (! skipSet.contains(genomeData.getGenomeId())) {
                RepGenome rep = new RepGenome(genomeData.getFid(), genomeData.getGenomeName(),
                        genomeData.getProtein());
                if (! genomeData.checkRepresented(sim) && ! repGenDb.checkSimilarity(rep, sim))
                    repGenDb.addRep(rep);
                batchCount--;
                if (batchCount == 0) {
                    double rate = genomeCount * 1000.0 / ((double) (System.currentTimeMillis() - start));
                    log.info("{} genomes out of {} processed for rep{}. {} genomes/second.",
                            genomeCount, this.genomeList.size(), sim, rate);
                    batchCount = this.batchSize;
                }
            }
        }
    }

    /**
     * Write the RepFinder file used to find close genomes.
     *
     * @throws FileNotFoundException
     */
    protected void writeRepFinder() throws FileNotFoundException {
        try (PrintWriter repFinderStream = new PrintWriter(new File(this.outDir, "repFinder.db"))) {
            // This requires two passes through the genome list.  The first lists the representing
            // genomes, the second lists the represented genomes.
            log.info("Creating repFinder.db using RepGen.{}.", this.lastRepGen.getThreshold());
            repFinderStream.println("genome_id\tgenome_name\tquality\trep_id\tgenetic_code\tseed_prot");
            // Loop through the genomes in the master list, writing the ones in the repgen set.
            for (ProteinData genome : this.genomeList) {
                if (this.lastRepGen.get(genome.getGenomeId()) != null) {
                    this.printRepFinderLine(repFinderStream, genome, genome.getGenomeId());
                }
            }
            // Write the spacer.
            repFinderStream.println("//");
            log.info("Representative genomes written.  Writing residual genomes.");
            // Loop through the saved lists of represented genomes, in statMap order.
            QualityCountMap<String> lastStatMap = this.statMaps.get(this.repGenSets.size() - 1);
            for (String groupId : lastStatMap.bestKeys()) {
                List<ProteinData> repFinderSet = this.repFinderSets.get(groupId);
                if (repFinderSet != null) {
                    for (ProteinData genome : this.repFinderSets.get(groupId)) {
                        this.printRepFinderLine(repFinderStream, genome, groupId);
                    }
                }
            }
        }
    }

    /**
     * Write the representative FASTA used by the PATRIC reference-genome computer.  This consists
     * of the seed proteins for the two highest-quality genomes in each repgen set at the highest
     * level.
     *
     * @throws IOException
     */
    protected void writeRefGenomeFasta() throws IOException {
        log.info("Writing FASTA file for PATRIC reference-genome computer.");
        // Get the highest-level repgen set.
        int repIdx = this.repGenSets.size() - 1;
        RepGenomeDb repGenSet = this.repGenSets.get(repIdx);
        // Create a count map for the repgen sets.  We output the first two in each.
        CountMap<String> repCounts = new CountMap<>();
        // Open the output file.
        try (FastaOutputStream fastaStream = new FastaOutputStream(new File(this.outDir, "refGenomes.fa"))) {
            int seqsOut = 0;
            int genomesChecked = 0;
            int totalGenomes = this.genomeList.size();
            // Loop through the protein data objects.  These are presented in quality order.
            log.info("Scanning protein database with {} genomes.", totalGenomes);
            for (ProteinData protein : this.genomeList) {
                // Get the representative.
                var representative = protein.getRepresentation(repGenSet);
                String repId = representative.getGenomeId();
                if (repCounts.count(repId) <= 2) {
                    // Here we have the first or second represented genome for this group.  Write it
                    // out.
                    Sequence seq = new Sequence(protein.getGenomeId(), protein.getGenomeName(),
                            protein.getProtein());
                    fastaStream.write(seq);
                    seqsOut++;
                }
                genomesChecked++;
                if (log.isInfoEnabled() && genomesChecked % 1000 == 0)
                    log.info("{} of {} genomes checked, {} output from {} groups.", genomesChecked,
                            totalGenomes, seqsOut, repCounts.size());
            }
        }
    }

    /**
     * Write the repgen protein FASTA files, the four-column tables, and the stat files.
     *
     * @throws IOException
     */
    protected void writeProteinFasta() throws IOException {
        for (int i = 0; i < this.repGenSets.size(); i++) {
            RepGenomeDb repGenSet = this.repGenSets.get(i);
            QualityCountMap<String> statMap = this.statMaps.get(i);
            String fileName = String.format("rep%d", repGenSet.getThreshold());
            log.info("Writing FASTA and statistics files for {}.", fileName);
            try (FastaOutputStream fastaStream = new FastaOutputStream(new File(this.outDir, fileName + ".faa"));
                    PrintWriter statsStream = new PrintWriter(new File(this.outDir, fileName + ".stats.tbl"));
                    PrintWriter seqStream = new PrintWriter(new File(this.outDir, fileName + ".seqs.tbl"))) {
                // Write the header for the stats file and the four-column table.
                statsStream.println("rep_id\trep_name\trating\tmembers");
                seqStream.println(FOUR_COLUMN_HEADER);
                // Loop through the groups from largest to smallest.
                for (String groupId : statMap.bestKeys()) {
                    RepGenome rep = repGenSet.get(groupId);
                    String groupName = rep.getName();
                    ProteinData genomeData = this.genomeList.getGenome(groupId);
                    fastaStream.write(new Sequence(groupId, groupName, rep.getProtein()));
                    statsStream.format("%s\t%s\t%s\t%d%n", groupId, groupName,
                            genomeData.getRating(), statMap.good(groupId));
                    this.writeFourCol(seqStream, genomeData);
                }
            }
        }
        // Write the master four-column table.
        try (PrintWriter fourColStream = new PrintWriter(new File(this.outDir, "all.seqs.tbl"))) {
            fourColStream.println(FOUR_COLUMN_HEADER);
            for (ProteinData genomeData : this.genomeList)
                this.writeFourCol(fourColStream, genomeData);
        }
    }

    /**
     * Write a record to a four-column table.
     *
     * @param seqStream		output stream
     * @param genomeData	genome descriptor for the genome to write
     */
    private void writeFourCol(PrintWriter seqStream, ProteinData genomeData) {
        seqStream.format("%s\t%s\t%s\t%s\t%s%n", genomeData.getGenomeId(),
                genomeData.getGenomeName(), genomeData.getProtein(),
                genomeData.getSsuSequence(), genomeData.getDna());
    }

    /**
     * Assign the genomes to their repgen sets.  Part of this is creating an output list file for
     * each set.  We also create count maps for the stat files we produce later.  For the last
     * repgen set, we save the assignments.
     *
     * @throws FileNotFoundException
     */
    protected void writeListFiles() throws FileNotFoundException {
        List<PrintWriter> listFiles = new ArrayList<>(this.repGenSets.size());
        this.statMaps = new ArrayList<>(this.repGenSets.size());
        // Here we save the groups for the last repgen set.
        this.lastRepGen = this.repGenSets.get(this.repGenSets.size() - 1);
        this.repFinderSets = new HashMap<>(this.lastRepGen.size());
        try {
            // First, create the output files and count maps and write the header lines.
            for (RepGenomeDb repGen : this.repGenSets) {
                File listFileName = new File(this.outDir, repGen.getListFileName());
                PrintWriter listFile = new PrintWriter(listFileName);
                listFile.println("genome_id\tgenome_name\tdomain\tgenus\tspecies\trep_id\tscore\tdistance");
                listFiles.add(listFile);
                this.statMaps.add(new QualityCountMap<>());
            }
            // These are used to generate progress messages.
            int genomeCount = 0;
            int batchCount = this.batchSize;
            long start = System.currentTimeMillis();
            for (ProteinData genome : this.genomeList) {
                // Get the descriptive part of the output line.
                String header = String.format("%s\t%s\t%s\t%s\t%s", genome.getGenomeId(),
                        genome.getGenomeName(), genome.getDomain(), genome.getGenus(),
                        genome.getSpecies());
                // Place the genome for each repgen set and output the appropriate line.
                for (int i = 0; i < this.repGenSets.size(); i++) {
                    // Write out the representation result.
                    RepGenomeDb repGenSet = this.repGenSets.get(i);
                    RepGenomeDb.Representation rep = genome.getRepresentation(repGenSet);
                    listFiles.get(i).format("%s\t%s\t%d\t%4.2f%n", header, rep.getGenomeId(),
                            rep.getSimilarity(), rep.getDistance());
                    // Update the count maps.
                    if (rep.isRepresented()) {
                        String repGenomeId = rep.getGenomeId();
                        this.statMaps.get(i).setGood(repGenomeId);
                        if (repGenSet == this.lastRepGen  && ! repGenomeId.contentEquals(genome.getGenomeId())) {
                            // Here we need to save the group info for the repfinder output.
                            this.repFinderSets.computeIfAbsent(repGenomeId, k -> new ArrayList<ProteinData>()).add(genome);
                        }
                        // Save the representation for the good-genome report.
                        genome.setRepresentation(repGenSet, rep);
                    } else {
                        this.statMaps.get(i).setBad(rep.getGenomeId());
                    }
                }
                genomeCount++;
                batchCount--;
                if (batchCount == 0) {
                    double rate = genomeCount * 1000.0 / ((double) (System.currentTimeMillis() - start));
                    log.info("{} genomes out of {} placed. {} genomes/second.",
                            genomeCount, this.genomeList.size(), rate);
                    batchCount = this.batchSize;
                }
            }
        } finally {
            // Close all the files.
            for (PrintWriter listFile : listFiles)
                listFile.close();
        }
    }

    /**
     * Write the seed-protein finder file.
     *
     * @throws IOException
     */
    protected void writeSeedProt() throws IOException {
        File seedFastaFile = new File(this.outDir, "seedProt.fa");
        log.info("Writing seed-finder FASTA to {}.", seedFastaFile);
        try (FastaOutputStream seedFasta = new FastaOutputStream(seedFastaFile)) {
            RepGenomeDb smallRepGen = this.repGenSets.get(0);
            for (RepGenome genome : smallRepGen) {
                ProteinData genomeDatum = this.genomeList.getGenome(genome.getGenomeId());
                Sequence seq = new Sequence(genomeDatum.getFid(), genomeDatum.getDomain(),
                        genomeDatum.getProtein());
                seedFasta.write(seq);
            }
        }
    }

    /**
     * Save all the RepGen sets.
     *
     * @throws IOException
     */
    protected void saveRepGenSets() throws IOException {
        for (RepGenomeDb repGen : this.repGenSets) {
            int score = repGen.getThreshold();
            String saveName = String.format("rep%d.ser", score);
            File saveFile = new File(this.outDir, saveName);
            log.info("Writing RepGen.{} to {}.  {} genomes in set.", score, saveFile,
                    repGen.size());
            repGen.save(saveFile);
        }
    }

    /**
     * Create the FASTA files for binning.
     *
     * @throws IOException
     */
    protected void createFastaFiles() throws IOException {
        try (FastaOutputStream dnaFasta = new FastaOutputStream(new File(this.outDir, "PhenTrnaSyntAlph.fa"));
                FastaOutputStream protFasta = new FastaOutputStream(new File(this.outDir, "allProts.fa"));
                FastaOutputStream ssuFasta = new FastaOutputStream(new File(this.outDir, "allSsu.fa"))) {
            log.info("Writing seed protein and SSU FASTA files.");
            for (ProteinData genomeData : this.genomeList) {
                String fid = genomeData.getFid();
                String gid = genomeData.getGenomeId();
                String name = genomeData.getGenomeName();
                // The DNA file expects the comment to be the ID, a tab, and the name.
                Sequence gSequence = new Sequence(fid, gid + "\t" + name, genomeData.getDna());
                dnaFasta.write(gSequence);
                // The protein file and SSU file just want a name.
                gSequence.setComment(name);
                gSequence.setSequence(genomeData.getProtein());
                protFasta.write(gSequence);
                gSequence.setSequence(genomeData.getSsuSequence());
                ssuFasta.write(gSequence);
            }
        }
    }

    /**
     * Validate the output directory and input file parameters.
     *
     * @throws IOException
     */
    protected void checkParms() throws IOException {
        File myOutDir = this.getOutDir();
        if (! myOutDir.exists()) {
            // Insure we have an output directory.
            log.info("Creating output directory {}.", myOutDir);
            if (! myOutDir.mkdir())
                throw new IOException("Could not create output directory " + myOutDir);
        } else if (! myOutDir.isDirectory()) {
            throw new FileNotFoundException("Invalid output directory " + myOutDir);
        }
    }
    /**
     * Read the input file and initialize the protein data genome list.  This also writes an output file
     * containing the genomes rejected for bad lineages.
     *
     * @param inFile	input file containing the evaluation data
     *
     * @return the protein data genome list
     *
     * @throws IOException
     * @throws UnsupportedEncodingException
     */
    protected ProteinDataFactory initializeProteinData(File inFile) throws IOException, UnsupportedEncodingException {
        // Initialize the genome list.
        this.genomeList = new ProteinDataFactory();
        // We use this to decide when to output progress messages.
        int batchCount = this.batchSize;
        // Read in the input file and get the protein data we need.
        try (TabbedLineReader inStream = new TabbedLineReader(inFile);
                PrintWriter saveStream = new PrintWriter(new File(this.outDir, "missing.taxons.tbl"))) {
            // Start the file where we save the genomes rejected for bad taxonomy data.
            saveStream.println("genome_id\tgenome_name\tscore\tlineage");
            // Loop through the file.
            for (TabbedLineReader.Line line : inStream) {
                // Accumulate statistics for this genome.
                this.genomeList.analyze(line);
                // Only process good genomes.
                if (line.getFlag(GenomeEval.GOOD_COL)) {
                    // Add the good genome to the list.
                    final String genomeId = line.get(GenomeEval.GENOME_COL);
                    final String genomeName = line.get(GenomeEval.NAME_COL);
                    final String lineageString = line.get(GenomeEval.LINEAGE_COL);
                    final double score = line.getDouble(GenomeEval.SCORE_COL);
                    boolean ok = this.genomeList.addGenome(genomeId, genomeName, lineageString, score);
                    if (! ok)
                        saveStream.format("%s\t%s\t%6.2f\t%s%n", genomeId, genomeName, score, lineageString);
                    // Show periodic progress.
                    batchCount--;
                    if (batchCount == 0) {
                        log.info("{} genomes found.", this.genomeList.size());
                        batchCount = this.batchSize;
                    }
                }
            }
        }
        log.info("{} good genomes collected.", this.genomeList.size());
        this.finishGenomeList(this.genomeList);
        // Remove the bad SSUs, if necessary.
        this.genomeList.prune(this.minRating);
        log.info("{} good genomes remaining after finishing.", this.genomeList.size());
        // Finally, the stats file.
        File statsFile = new File(this.outDir, "bvbrc.final.stats.tbl");
        this.genomeList.writeStats(statsFile);
        return this.genomeList;
    }

    /**
     * Finish the genome list.  This involves analyzing the SSUs and filling in the protein sequences.
     *
     * @throws UnsupportedEncodingException
     * @throws IOException
     */
    protected abstract void finishGenomeList(ProteinDataFactory gList)
            throws UnsupportedEncodingException, IOException;

    /**
     * Write a single line of the repFinder.db file.
     *
     * @param repFinderStream	output stream
     * @param genome			proteinData object for the genome to write
     * @param repId				associated representative genome ID
     */
    private void printRepFinderLine(PrintWriter repFinderStream, ProteinData genome, String repId) {
        repFinderStream.format("%s\t%s\t%4.4f\t%s\t%d\t%s%n", genome.getGenomeId(),
                genome.getGenomeName(), genome.getScore(), repId, genome.getGeneticCode(),
                genome.getProtein());
    }

    /**
     * Create the array of RepGen databases.
     *
     * @param size	estimated number needed
     */
    @Override
    public void initRepGenSets(int size) {
        this.repGenSets = new ArrayList<>(size);
    }

    /**
     * Store a new RepGen database in the list.
     *
     * @param repGenomeDb	RepGen database to store
     */
    @Override
    public void addRepGenSet(RepGenomeDb repGenomeDb) {
        this.repGenSets.add(repGenomeDb);
    }

    /**
     * Write the final list of good genomes.
     *
     * @throws IOException
     */
    protected void writeGoodGenomes() throws IOException {
        // Build the headers for the representatives.
        int[] repLevels = this.repGenSets.stream().mapToInt(x -> x.getThreshold()).toArray();
        String headers = Arrays.stream(repLevels).mapToObj(x -> String.format("rep%d", x))
                .collect(Collectors.joining("\t"));
        File goodFile = new File(this.outDir, "bvbrc.good.tbl");
        try (PrintWriter writer = new PrintWriter(goodFile)) {
            // Write the header.
            writer.println("genome_id\tname\tdomain\tgenus\tspecies\tscore\trating\t" + headers);
            // We'll build the output in here.
            TextStringBuilder buffer = new TextStringBuilder(150);
            for (ProteinData genome : this.genomeList) {
                // Write the static data.
                buffer.append(genome.getGenomeId()).append('\t');
                buffer.append(genome.getGenomeName()).append('\t');
                buffer.append(genome.getDomain()).append('\t');
                buffer.append(genome.getGenus()).append('\t');
                buffer.append(genome.getSpecies()).append('\t');
                buffer.append("%8.4f\t", genome.getScore());
                buffer.append(genome.getRating().toString());
                // Loop through the repgens.  Here we append the tab at the front.
                for (int repLevel : repLevels)
                    buffer.append('\t').append(genome.getRepGenome(repLevel));
                // Write the line.
                writer.println(buffer.toString());
                buffer.clear();
            }
            writer.flush();
        }
    }

    /**
     * Build new GTO directories for the repgen sets.
     *
     * @throws IOException
     */
    public void setupGTOs() throws IOException {
        for (RepGenomeDb repdb : this.repGenSets) {
            File repDir = computeRepGtoDir(repdb);
            if (! repDir.isDirectory()) {
                log.info("Creating GTO output directory {}.", repDir);
                FileUtils.forceMkdir(repDir);
            } else {
                log.info("Clearing GTO output directory {}.", repDir);
                FileUtils.cleanDirectory(repDir);
            }
        }
    }

    /**
     * @return the name of the GTO output directory for a repgen set
     *
     * @param repdb		repgen set database
     */
    private File computeRepGtoDir(RepGenomeDb repdb) {
        return new File(this.outDir, "GTO" + Integer.toString(repdb.getThreshold()));
    }

    /**
     * Compute where all of the representative genomes go in the output directories.
     */
    public void placeGTOs() throws IOException {
        int outCount = 0;
        for (RepGenomeDb repdb : this.repGenSets) {
            // Loop through this repgen set, connecting the genome IDs to this directory.
            for (var rep : repdb) {
                File repDir = this.computeRepGtoDir(repdb);
                String genomeId = rep.getGenomeId();
                this.queueGtoDownload(repDir, genomeId);
                outCount++;
            }
        }
        log.info("{} distinct genomes will be downloaded to {} GTO files.", repHash.size(), outCount);

    }

    /**
     * Copy GTOs from old repgen GTO directories to new ones.
     *
     * @param inDir		input directory containing the old GTOs
     *
     * @throws IOException
     */
    public void copyGTOs(File inDir) throws IOException {
        // We count the number of genomes copied and the number of additional copies required.
        int outCount = 0;
        int copyCount = 0;
        // This map will tell us where to find each genome that is already downloaded.
        Map<String, File> oldMap = new HashMap<>();
        // Get the GTO file names.
        var gtoFiles = FileUtils.listFiles(inDir, GTO_FILE_FILTER, GTO_DIR_FILTER);
        log.info("{} existing representative GTO files found in {}.", gtoFiles.size(), inDir);
        // Now we must loop through all the GTOs and put them in the map.
        for (File gtoFile : gtoFiles) {
            // The genome ID is the prefix of the file name.  The file filter pattern insures
            // that all we need to do is chop off the ".gto" at the end.
            String genomeId = StringUtils.removeEnd(gtoFile.getName(), ".gto");
            if (! oldMap.containsKey(genomeId)) {
                // Here the genome is a new one, so we save its location.
                oldMap.put(genomeId, gtoFile);
            }
        }
        log.info("{} unique representative genomes found.", oldMap.size());
        // Now we process the repgen sets.  For each, we copy the genomes that already exist,
        // and queue the others for downloading.
        for (RepGenomeDb repdb : this.repGenSets) {
            // Loop through this repgen set.
            for (var rep : repdb) {
                File repDir = this.computeRepGtoDir(repdb);
                String genomeId = rep.getGenomeId();
                if (oldMap.containsKey(genomeId)) {
                    // Here we have to copy.
                    File inFile = oldMap.get(genomeId);
                    File outFile = new File(repDir, genomeId + ".gto");
                    log.info("Copying {} ({}) from {} to {}.", genomeId, rep.getName(), inFile, outFile);
                    FileUtils.copyFile(inFile, outFile);
                    copyCount++;
                } else {
                    // Here we have to queue for download.
                    this.queueGtoDownload(repDir, genomeId);
                    outCount++;
                }
            }
        }
        log.info("{} genomes copied, {} queued for download to {} files.", copyCount, this.repHash.size(), outCount);
    }

    /**
     * Specify that a genome needs to be downloaded to the specified directory.  This insures there
     * is a repHash entry for the genome and adds the directory to its directory list.
     *
     * @param repDir	target output directory
     * @param genomeId	genome to download
     */
    private void queueGtoDownload(File repDir, String genomeId) {
        List<File> dirList = this.repHash.computeIfAbsent(genomeId,
                x -> new ArrayList<File>(this.repGenSets.size()));
        dirList.add(repDir);
    }


    /**
     * Save the new repgen GTOs to directories.
     *
     * @throws IOException
     */
    public void saveNewGTOs() throws IOException {
        // Connect to PATRIC.
        P3CursorConnection p3 = new P3CursorConnection();
        // Loop through all the genomes doing the downloads.
        int gCount = 0;
        for (var repEntry : repHash.entrySet()) {
            String genomeId = repEntry.getKey();
            gCount++;
            log.info("Downloading {} of {}: {}.", gCount, repHash.size(), genomeId);
            Genome genome = P3Genome.load(p3, genomeId, P3Genome.Details.FULL);
            String gName = genomeId + ".gto";
            for (File myOutDir : repEntry.getValue())
                genome.save(new File(myOutDir, gName));
        }
    }

    /**
     * @return the number of genomes in this run
     */
    public int genomeCount() {
        return this.genomeList.size();
    }

    /**
     * @return TRUE if the specified genome is present in this run
     *
     * @param genomeId 	ID of the genome in question
     */
    public boolean isPresent(String genomeId) {
        return this.genomeList.contains(genomeId);
    }

    /**
     * @return TRUE if GTOs should be downloaded
     */
    public boolean isGTOsRequested() {
        return (this.gtoMode == GtoScheme.YES || this.gtoMode == GtoScheme.ONLY);
    }

    /**
     * @return TRUE if the repgen sets should be rebuilt
     */
    public boolean isFullRegenRequired() {
        return (this.gtoMode == GtoScheme.YES || this.gtoMode == GtoScheme.NO);
    }

    /**
     * Save the list of repgen set levels in the protein data factory.
     */
    public void saveRepLevels() {
        int[] levels = this.repGenSets.stream().mapToInt(x -> x.getThreshold()).toArray();
        this.genomeList.setRepLevels(levels);
    }

    /**
     * Write out the genome report from the protein data factory.
     *
     * @throws IOException
     */
    public void writeGenomeReport() throws IOException {
        File gReportFile = new File(this.outDir, "bvbrc.genome.tbl");
        this.genomeList.writeReport(gReportFile);
    }

}
