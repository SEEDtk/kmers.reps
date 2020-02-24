/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
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
import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ICommand;
import org.theseed.utils.IntegerList;

/**
 * This command processes the results from a quality run.  It expects as input the sorted results
 * file containing the lineage, the quality data, and the good-genome flag.  It will create
 * the master seed protein FASTA file, the master binning BLAST database, and the representative-genome
 * databases at level 10, 50, 100, and 200.
 *
 * The positional parameters are the name of the output directory and the name of the input file.
 *
 * The following files are put into the output directory.
 *
 * The following command-line options are supported.
 *
 * -v	show more detailed log messages
 * -b	batch size for PATRIC queries
 *
 * --repSizes	comma-delimited list of score limits for each repGen set
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor extends BaseProcessor implements ICommand {

    // FIELDS
    /** manager for key data about each genome and its seed protein */
    private ProteinDataFactory genomeList;
    /** list of RepGen set thresholds */
    private IntegerList repSizes;

    // COMMAND-LINE OPTIONS

    /** number of genomes per batch when retrieving sequences */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "200",
            usage = "number of genomes per batch")
    private int batchSize;

    /** list of RepGen set scores */
    @Option(name = "--repSizes", aliases = { "--repScores", "--repSims" }, metaVar = "50,100,200,250",
            usage = "comma-delimited list of minimum scores for each RepGen set to create")
    private void setRepSizes(String repSizeString) {
        this.repSizes = new IntegerList(repSizeString);
    }

    /** output directory */
    @Argument(index = 0, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    /** input tab-delimited file of genomes */
    @Argument(index = 1, metaVar = "inFile.tbl", usage = "input file", required = true)
    private File inFile;

    @Override
    protected void setDefaults() {
        this.batchSize = 500;
        this.repSizes = new IntegerList(10, 50, 100, 200);
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.outDir.exists()) {
            // Insure we have an output directory.
            log.info("Creating output directory {}.", this.outDir);
            if (! this.outDir.mkdir())
                throw new IOException("Could not create output directory " + this.outDir);
        } else if (! this.outDir.isDirectory()) {
            throw new FileNotFoundException("Invalid output directory " + this.outDir);
        }
        if (! this.inFile.canRead())
            throw new FileNotFoundException(this.inFile + " is not found or unreadable.");
        return true;
    }

    @Override
    public void run() {
        try {
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
            // We need to create the FASTA files for the seed protein list and the
            // binning BLAST database.  We do that here.
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
            // Create a list of RepGenomeDb objects, one for each score limit.
            log.info("Creating repgen sets.");
            List<RepGenomeDb> repGenSets = new ArrayList<RepGenomeDb>(this.repSizes.size());
            for (int repSize : this.repSizes) {
                repGenSets.add(new RepGenomeDb(repSize, ProteinDataFactory.SEED_FUNCTION));
            }
            log.info("Sorting genomes into repgen sets.");
            int genomeCount = 0;
            long start = System.currentTimeMillis();
            batchCount = this.batchSize;
            for (ProteinData genomeData : this.genomeList) {
                RepGenome rep = new RepGenome(genomeData.getFid(), genomeData.getGenomeName(),
                        genomeData.getProtein());
                for (RepGenomeDb repGen : repGenSets) {
                    if (! repGen.checkSimilarity(rep, repGen.getThreshold())) {
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
            // Save all the repgen sets.
            for (RepGenomeDb repGen : repGenSets) {
                int score = repGen.getThreshold();
                String saveName = String.format("rep%d.ser", score);
                File saveFile = new File(this.outDir, saveName);
                log.info("Writing RepGen.{} to {}.  {} genomes in set.", score, saveFile,
                        repGen.size());
                repGen.save(saveFile);
            }
            // Write out the protein Fasta file for the first set.  This is used to find
            // seed proteins.
            File seedFastaFile = new File(this.outDir, "seedProt.fa");
            log.info("Writing seed-finder FASTA to {}.", seedFastaFile);
            try (FastaOutputStream seedFasta = new FastaOutputStream(seedFastaFile)) {
                RepGenomeDb smallRepGen = repGenSets.get(0);
                for (RepGenome genome : smallRepGen) {
                    ProteinData genomeDatum = this.genomeList.getGenome(genome.getGenomeId());
                    Sequence seq = new Sequence(genomeDatum.getFid(), genomeDatum.getDomain(),
                            genomeDatum.getProtein());
                    seedFasta.write(seq);
                }
            }
            // Now we need to assign the genomes to their repgen sets.  Part of this is creating
            // an output text file for each set.  We also create count maps for the stat files
            // we produce later.  For the last repgen set, we save the assignments.
            log.info("Assigning genomes to repgen sets.");
            List<PrintWriter> listFiles = new ArrayList<PrintWriter>(repGenSets.size());
            List<QualityCountMap<String>> statMaps = new ArrayList<QualityCountMap<String>>(repGenSets.size());
            // Here we save the groups for the last repgen set.
            RepGenomeDb lastRepGen = repGenSets.get(repGenSets.size() - 1);
            Map<String, List<ProteinData>> repFinderSets = new HashMap<String, List<ProteinData>>(lastRepGen.size());
            try {
                // First, create the output files and count maps and write the header lines.
                for (RepGenomeDb repGen : repGenSets) {
                    File listFileName = new File(this.outDir,
                            String.format("rep%d.list.tbl", repGen.getThreshold()));
                    PrintWriter listFile = new PrintWriter(listFileName);
                    listFile.println("genome_id\tgenome_name\tdomain\tgenus\tspecies\trep_id\tscore\tdistance");
                    listFiles.add(listFile);
                    statMaps.add(new QualityCountMap<String>());
                }
                // These are used to generate progress messages.
                genomeCount = 0;
                batchCount = this.batchSize;
                start = System.currentTimeMillis();
                for (ProteinData genome : this.genomeList) {
                    // Get the kmers for this genome.
                    ProteinKmers protein = new ProteinKmers(genome.getProtein());
                    // Get the descriptive part of the output line.
                    String header = String.format("%s\t%s\t%s\t%s\t%s", genome.getGenomeId(),
                            genome.getGenomeName(), genome.getDomain(), genome.getGenus(),
                            genome.getSpecies());
                    // Place the genome for each repgen set and output the appropriate line.
                    for (int i = 0; i < repGenSets.size(); i++) {
                        // Write out the representation result.
                        RepGenomeDb repGenSet = repGenSets.get(i);
                        RepGenomeDb.Representation rep = repGenSet.findClosest(protein);
                        listFiles.get(i).format("%s\t%s\t%d\t%4.2f%n", header, rep.getGenomeId(),
                                rep.getSimilarity(), rep.getDistance());
                        // Update the count maps.
                        if (rep.isRepresented()) {
                            String repGenomeId = rep.getGenomeId();
                            statMaps.get(i).setGood(repGenomeId);
                            if (repGenSet == lastRepGen  && ! repGenomeId.contentEquals(genome.getGenomeId())) {
                                // Here we need to save the group info for the repfinder output.
                                repFinderSets.computeIfAbsent(repGenomeId, k -> new ArrayList<ProteinData>()).add(genome);
                            }
                        } else {
                            statMaps.get(i).setBad(rep.getGenomeId());
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
            // Now we write the protein FASTA files and the stats files.
            for (int i = 0; i < repGenSets.size(); i++) {
                RepGenomeDb repGenSet = repGenSets.get(i);
                QualityCountMap<String> statMap = statMaps.get(i);
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
            // Now we produce the repFinder file used to find close genomes.  This requires two passes
            // through the genome list.  The first lists the representing genomes, the second lists
            // the represented genomes.
            try (PrintWriter repFinderStream = new PrintWriter(new File(this.outDir, "repFinder.db"))) {
                log.info("Creating repFinder.db using RepGen.{}.", lastRepGen.getThreshold());
                repFinderStream.println("genome_id\tgenome_name\tquality\trep_id\tgenetic_code\tseed_prot");
                // Loop through the genomes in the master list, writing the ones in the repgen set.
                for (ProteinData genome : this.genomeList) {
                    if (lastRepGen.get(genome.getGenomeId()) != null) {
                        this.printRepFinderLine(repFinderStream, genome, genome.getGenomeId());
                    }
                }
                // Write the spacer.
                repFinderStream.println("//");
                log.info("Representative genomes written.  Writing residual genomes.");
                // Loop through the saved lists of represented genomes, in statMap order.
                QualityCountMap<String> lastStatMap = statMaps.get(repGenSets.size() - 1);
                for (String groupId : lastStatMap.bestKeys()) {
                    List<ProteinData> repFinderSet = repFinderSets.get(groupId);
                    if (repFinderSet != null) {
                        for (ProteinData genome : repFinderSets.get(groupId)) {
                            this.printRepFinderLine(repFinderStream, genome, groupId);
                        }
                    }
                }
            }
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void printRepFinderLine(PrintWriter repFinderStream, ProteinData genome, String repId) {
        repFinderStream.format("%s\t%s\t%4.4f\t%s\t%d\t%s%n", genome.getGenomeId(),
                genome.getGenomeName(), genome.getScore(), repId, genome.getGeneticCode(),
                genome.getProtein());
    }

}
