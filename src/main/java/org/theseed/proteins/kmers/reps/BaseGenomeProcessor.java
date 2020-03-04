package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.counters.QualityCountMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.ProteinData;
import org.theseed.proteins.ProteinDataFactory;
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

    /** number of genomes per batch when retrieving sequences */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "200", usage = "number of genomes per batch")
    private int batchSize;
    /** output directory */
    @Argument(index = 0, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;
    /** input tab-delimited file of genomes */
    @Argument(index = 1, metaVar = "inFile.tbl", usage = "input file", required = true)
    private File inFile;


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

    /**
     * @return the input file
     */
    protected File getInFile() {
        return inFile;
    }

    /**
     * Update the batch size.
     *
     * @param batchSize the proposed new batch size
     */
    protected void setBatchSize(int batchSize) {
        this.batchSize = batchSize;
    }

    public BaseGenomeProcessor() {
        super();
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
     * Write the protein FASTA files and the stat files.
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
                    PrintWriter statsStream = new PrintWriter(new File(this.outDir, fileName + ".stats.tbl"))) {
                // Write the header for the stats file.
                statsStream.println("rep_id\trep_name\tmembers\toutliers");
                // Loop through the groups from largest to smallest.
                for (String groupId : statMap.bestKeys()) {
                    RepGenome rep = repGenSet.get(groupId);
                    String groupName = rep.getName();
                    fastaStream.write(new Sequence(groupId, groupName, rep.getProtein()));
                    statsStream.format("%s\t%s\t%d\t%d%n", groupId, groupName, statMap.good(groupId),
                            statMap.bad(groupId));
                }
            }
        }
    }

    /**
     * Assign the genomes to their repgen sets.  Part of this is creating an output test file for
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
                FastaOutputStream protFasta = new FastaOutputStream(new File(this.outDir, "allProts.fa"))) {
            log.info("Writing seed protein FASTA files.");
            for (ProteinData genomeData : this.genomeList) {
                String fid = genomeData.getFid();
                String gid = genomeData.getGenomeId();
                String name = genomeData.getGenomeName();
                // The DNA file expects the comment to be the ID, a tab, and the name.
                Sequence gSequence = new Sequence(fid, gid + "\t" + name, genomeData.getDna());
                dnaFasta.write(gSequence);
                // The protein file wants just a name.
                gSequence.setComment(name);
                gSequence.setSequence(genomeData.getProtein());
                protFasta.write(gSequence);
            }
        }
    }

    /**
     * Read the input file and initialize the protein data genome list.
     *
     * @throws IOException
     * @throws UnsupportedEncodingException
     */
    protected void initializeProteinData() throws IOException, UnsupportedEncodingException {
        // Initialize the genome list.
        this.genomeList = new ProteinDataFactory();
        // We use this to decide when to output progress messages.
        int batchCount = this.batchSize;
        // Read in the input file and get the protein data we need.
        try (TabbedLineReader inStream = new TabbedLineReader(this.inFile)) {
            // Compute the field indices.
            int idCol = inStream.findField("genome_id");
            int nameCol = inStream.findField("genome_name");
            int lineageCol = inStream.findField("taxon_lineage_ids");
            int scoreCol = inStream.findField("score");
            int goodCol = inStream.findField("Good Genome");
            // Loop through the file.
            for (TabbedLineReader.Line line : inStream) {
                // Only process good genomes.
                if (line.getFlag(goodCol)) {
                    this.genomeList.addGenome(line.get(idCol), line.get(nameCol),
                            line.get(lineageCol), line.getDouble(scoreCol));
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
        this.genomeList.finishList(this.batchSize);
        log.info("{} good genomes remaining after seed protein scan.", this.genomeList.size());
    }

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

}
