/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;

/**
 * This command will produce a report on the N closest representative genomes to each genome in a genome source.  The
 * default is to show the top 5.
 *
 * The output report will be produced on the standard output and will contain the similarity score and distance for
 * each input/representative genome pair along with the ID and name of each genome in the pair.
 *
 * The positional parameters should be the name of the representative-genome database and the name of the genome source
 * directory or file.
 *
 * The command-line options are as follows:
 *
 * -h	display the command-line usage
 * -v	display more frequent log messages
 * -o	output file for the report (if not STDOUT)
 * -n	number of close representatives to show (default 5)
 *
 * --min		minimum acceptable similarity to display (default 1)
 * --source		type of genome source (default DIR)
 * --para		use parallel processing to speed the search (default FALSE)
 *
 * @author Bruce Parrello
 *
 */
public class NeighborProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(NeighborProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** representative-genome database */
    private RepGenomeDb repDb;
    /** array of representative genomes */
    private RepGenome[] reps;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of input genome source")
    private GenomeSource.Type sourceType;

    /** number of close representatives to return */
    @Option(name = "--num", aliases = { "-n" }, metaVar = "2", usage = "number of close representatives to return per input genome")
    private int num;

    /** if specified, input genomes will be processed in parallel */
    @Option(name = "--para", usage = "if specified, parallel processing will be used to improve performance")
    private boolean parallelFlag;

    /** minimum acceptable similarity to display */
    @Option(name = "--min", metaVar = "100", usage = "minimum acceptable similarity to display")
    private int minSim;

    /** representative-genome database file */
    @Argument(index = 0, metaVar = "repXX.set", usage = "name of the representative-genome database file", required = true)
    private File repDbFile;

    /** genome source file or directory */
    @Argument(index = 1, metaVar = "inDir", usage = "name of the input genome source directory or file", required = true)
    private File inDir;

    /**
     * Utility object for data needed from the genome.
     */
    private class GenomeData {

        /** repgenome object for the seed protein */
        private RepGenome seedProtein;
        /** genome size */
        private int dnaSize;

        protected GenomeData(Genome genome) {
            this.seedProtein = NeighborProcessor.this.repDb.getSeedProtein(genome);
            this.dnaSize = genome.getLength();
        }

    }

    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.num = 5;
        this.parallelFlag = false;
        this.minSim = 1;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify the close-genome count.
        if (this.num <= 0)
            throw new ParseFailureException("Close-genome count (num) must be at least 1.");
        // Verify the minimum similarity.
        if (this.minSim < 0)
            throw new ParseFailureException("Minimum similarity cannot be zero.");
        // Verify the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " does not exist.");
        // Verify that we can read the repgen database.
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Repgen database file " + this.repDbFile + " is not found or unreadable.");
        // Connect to the source.
        log.info("Connecting to genome source {} of type {}.", this.inDir, this.sourceType.toString());
        this.genomes = this.sourceType.create(this.inDir);
        // Load the repgen database.
        log.info("Loading repgen database from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Output the report header.
        writer.println("genome_id\tgenome_name\tsize\trep_id\trep_name\tsim\tdistance\tclosest\tsignal");
        // Get the list of incoming genome IDs.
        var genomeIDs = this.genomes.getIDs();
        log.info("{} genomes found in source.", genomeIDs.size());
        // Get an array of representative-genome objects.
        this.reps = this.repDb.all();
        // Loop through them, searching for representatives.
        Stream<String> idStream = genomeIDs.stream();
        if (this.parallelFlag)
            idStream = idStream.parallel();
        idStream.forEach(x -> this.processGenome(writer, x));
    }

    /**
     * This method computes the N closest represenatives for a single genome.  We use a sorted list, and add each
     * new representation into the list, lopping off any at the end.  Then we write the results.
     *
     * @param writer	output report writer
     * @param genomeID	ID of input genome to process
     */
    private void processGenome(PrintWriter writer, String genomeID) {
        // Get the incoming genome's representation object and its size.
        GenomeData gData = this.getGenomeData(genomeID);
        RepGenome seedProtein = gData.seedProtein;
        if (seedProtein == null)
            log.error("No seed protein found for input genome {}.", genomeID);
        else {
            // Create the output queue.
            var found = new ArrayList<RepGenomeDb.Representation>(this.num + 1);
            // Loop through the representative genomes, keeping the best.
            long timer = System.currentTimeMillis();
            int repCount = 0;
            for (RepGenome rep : this.reps) {
                RepGenomeDb.Representation result = this.repDb.new Representation(rep, seedProtein);
                // Only keep it if the similarity is greater than the minimum
                final int sim = result.getSimilarity();
                if (sim >= this.minSim) {
                    // Add it to the found-results list.
                    int idx = 0;
                    final int n = found.size();
                    while (idx < n && found.get(idx).getSimilarity() > sim)
                        idx++;
                    // Now idx points to the place where we want to insert this representative.
                    if (idx < n || n < this.num) {
                        // Here this new result goes in the queue.  Either it is being added to the middle, or the
                        // queue has space.
                        found.add(idx, result);
                        // Now, if the queue is full, see if we want to get rid of the last one.
                        if (found.size() > this.num) {
                            // Compare the reps at the break point.  We know there is a rep at the break point and at
                            // least one past it.
                            if (found.get(this.num-1).getSimilarity() > found.get(this.num).getSimilarity())
                                // Everything past the break point is lower, so delete them all.
                                while (found.size() > this.num)
                                    found.remove(this.num);
                        }
                    }
                }
                repCount++;
                if (log.isInfoEnabled() && System.currentTimeMillis() - timer >= 10000) {
                    timer = System.currentTimeMillis();
                    log.info("{} representatives processed, {} kept.", repCount, found.size());
                }
            }
            if (found.size() <= 0)
                log.warn("No close genomes found for {}.", genomeID);
            this.writeResults(writer, gData, found);
        }
    }

    /**
     * Write the close genomes for the specified input genome to the output report.
     *
     * @param writer		output report writer
     * @param gData			input genome descriptor
     * @param found			representative genome object
     */
    private synchronized void writeResults(PrintWriter writer, GenomeData gData, List<RepGenomeDb.Representation> found) {
        var seedProtein = gData.seedProtein;
        String prefix = seedProtein.getGenomeId() + "\t" + seedProtein.getName() + "\t" + Integer.toString(gData.dnaSize);
        if (found.size() <= 0)
            writer.println(prefix + "\t\t\t\t\t\t");
        else {
            String closeFlag = "Y";
            for (RepGenomeDb.Representation rep : found) {
                int sim = rep.getSimilarity();
                double dist = rep.getDistance();
                String repId = rep.getGenomeId();
                double signal = gData.dnaSize * (1.0 - dist);
                writer.println(prefix + "\t" + repId + "\t" + this.repDb.get(repId).getName() + "\t" +
                Integer.toString(sim) + "\t" + Double.toString(dist) + "\t" + closeFlag + "\t" + Double.toString(signal));
                closeFlag = " ";
            }
        }
        writer.flush();
    }

    /**
     * @return the data we need from the specified input genome
     *
     * @param genomeID	ID of the desired input genome
     */
    private GenomeData getGenomeData(String genomeID) {
        Genome genome = this.genomes.getGenome(genomeID);
        log.info("Processing genome {}.", genome);
        GenomeData retVal = this.new GenomeData(genome);
        return retVal;
    }

}
