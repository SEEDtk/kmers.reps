package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.text.TextStringBuilder;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.counters.CountMap;
import org.theseed.counters.GenomeEval;
import org.theseed.counters.QualityCountMap;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseProcessor;

/**
 * This is the base class for both GenomeProcessor and UpdateProcessor.  It performs the common tasks
 * in maintaining the various files that depend on genome evaluations.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseGenomeProcessor extends BaseProcessor implements IRepGenContainer {

    // FIELDS
    /** header line for four-column table containing SSUs */
    private static final String FOUR_COLUMN_HEADER = "genome_id\tgenome_name\tseed_protein\tssu_rna\tseed_dna";

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

    // COMMAND-LINE OPTIONS

    /** specifies whether bad-SSU genomes will be considered good */
    @Option(name = "--minLevel", usage = "minimum rating level to keep")
    private ProteinData.Rating minRating;

    /** number of genomes per batch when retrieving sequences */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "200", usage = "number of genomes per batch")
    private int batchSize;

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
    }

    /**
     * Set the base-class option defaults.
     */
    protected void setBaseDefaults() {
        this.batchSize = 500;
        this.minRating = ProteinData.Rating.SINGLE_SSU;
    }

    /**
     * Separate the genomes into RepGen sets.
     */
    protected void collateGenomes() {
        log.info("Sorting genomes into repgen sets.");
        int genomeCount = 0;
        long start = System.currentTimeMillis();
        int batchCount = this.batchSize;
        for (ProteinData genomeData : this.genomeList) {
            RepGenome rep = new RepGenome(genomeData.getFid(), genomeData.getGenomeName(),
                    genomeData.getProtein());
            for (RepGenomeDb repGen : this.repGenSets) {
                int sim = repGen.getThreshold();
                if (! genomeData.checkRepresented(sim) && ! repGen.checkSimilarity(rep, sim)) {
                    repGen.addRep(rep);
                }
            }
            batchCount--;
            genomeCount++;
            if (batchCount == 0) {
                double rate = genomeCount * 1000.0 / ((double) (System.currentTimeMillis() - start));
                log.info("{} genomes out of {} processed. {} genomes/second.",
                        genomeCount, this.genomeList.size(), rate);
                batchCount = this.batchSize;
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
        CountMap<String> repCounts = new CountMap<String>();
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
        List<PrintWriter> listFiles = new ArrayList<PrintWriter>(this.repGenSets.size());
        this.statMaps = new ArrayList<QualityCountMap<String>>(this.repGenSets.size());
        // Here we save the groups for the last repgen set.
        this.lastRepGen = this.repGenSets.get(this.repGenSets.size() - 1);
        this.repFinderSets = new HashMap<String, List<ProteinData>>(this.lastRepGen.size());
        try {
            // First, create the output files and count maps and write the header lines.
            for (RepGenomeDb repGen : this.repGenSets) {
                File listFileName = new File(this.outDir,
                        String.format("rep%d.list.tbl", repGen.getThreshold()));
                PrintWriter listFile = new PrintWriter(listFileName);
                listFile.println("genome_id\tgenome_name\tdomain\tgenus\tspecies\trep_id\tscore\tdistance");
                listFiles.add(listFile);
                this.statMaps.add(new QualityCountMap<String>());
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
        File outDir = this.getOutDir();
        if (! outDir.exists()) {
            // Insure we have an output directory.
            log.info("Creating output directory {}.", outDir);
            if (! outDir.mkdir())
                throw new IOException("Could not create output directory " + outDir);
        } else if (! outDir.isDirectory()) {
            throw new FileNotFoundException("Invalid output directory " + outDir);
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
        File statsFile = new File(this.outDir, "patric.stats.tbl");
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
        this.repGenSets = new ArrayList<RepGenomeDb>(size);
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
        File goodFile = new File(this.outDir, "patric.good.tbl");
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
     * Save the repgen GTOs to directories.
     *
     * @throws IOException
     */
    public void saveGTOs() throws IOException {
        // Create an output directory for each repgen set.  We will connect each representative genome to
        // all the directories in which it belongs.
        var repHash = new TreeMap<String, List<File>>();
        int outCount = 0;
        for (RepGenomeDb repdb : this.repGenSets) {
            File repDir = new File(this.outDir, "GTO" + Integer.toString(repdb.getThreshold()));
            if (! repDir.isDirectory()) {
                log.info("Creating GTO output directory {}.", repDir);
                FileUtils.forceMkdir(repDir);
            } else {
                log.info("Clearing GTO output directory {}.", repDir);
                FileUtils.cleanDirectory(repDir);
            }
            // Loop through this repgen set, connecting the genome IDs to this directory.
            for (var rep : repdb) {
                String genomeId = rep.getGenomeId();
                List<File> dirList = repHash.computeIfAbsent(genomeId, x -> new ArrayList<File>(this.repGenSets.size()));
                dirList.add(repDir);
                outCount++;
            }
        }
        log.info("{} distinct genomes will be downloaded to {} GTO files.", repHash.size(), outCount);
        // Connect to PATRIC.
        P3Connection p3 = new P3Connection();
        // Loop through all the genomes doing the downloads.
        int gCount = 0;
        for (var repEntry : repHash.entrySet()) {
            String genomeId = repEntry.getKey();
            gCount++;
            log.info("Downloading {} of {}: {}.", gCount, repHash.size(), genomeId);
            Genome genome = P3Genome.load(p3, genomeId, P3Genome.Details.FULL);
            String gName = genomeId + ".gto";
            for (File outDir : repEntry.getValue())
                genome.save(new File(outDir, gName));
        }
    }

}
