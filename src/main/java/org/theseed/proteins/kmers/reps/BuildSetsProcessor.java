package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This command will create neighbor sets of a specified size from a genome representation list file. The
 * file used is one of the repXXX.list.tbl files created during the repgen build process. The neighbor sets
 * are created by finding, for each representing genome, its N closest neighbors, where N is the desired
 * number of genomes to be output in each generated list. The generated lists will be written as separate
 * files in the output directory.
 * 
 * In this case, closeness is determined by the kmer score, in the "score" column of the input file. For
 * each genome, we read in its ID and name (genome_id and genome_name), its representative ID (rep_id),
 * and its score (score). Each genome is stored in a map of lists keyed on representative ID. After all
 * genomes are read, we output the first N genomes from each list. If a list has fewer than N genomes,
 * it is skipped.
 * 
 * The goal is to produce diverse sets of similar genomes for testing purposes. The input list file is
 * read from the standard input. There are no positional parameters.
 * 
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -D	output directory name (default "Sets" in the current directory)
 * -i   input file name (if not STDIN)
 *
 * --min    minimum number of genomes per set (default 1500)
 * --max    maximum number of genomes per set (default 2000)
 * --clear	erase the output directory before processing
 * 
 * @author Bruce Parrello
 */

public class BuildSetsProcessor extends BaseMultiReportProcessor {

    // FIELDS
    /** logging facility */ 
    private static final Logger log = LoggerFactory.getLogger(BuildSetsProcessor.class);
    /** map of representative IDs to lists of genomes */
    private Map<String, List<GenomeRecord>> repMap;
    /** index of the input genome ID column */
    private int idCol;
    /** index of the input genome name column */
    private int nameCol;
    /** index of the input representative ID column */
    private int repCol;
    /** index of the input score column */
    private int scoreCol;
    /** input file stream */
    private TabbedLineReader inStream;

    // COMMAND-LINE OPTIONS

    /** minimum number of genomes per set */
    @Option(name = "--min", metaVar = "2000", usage = "minimum number of genomes per output set")
    private int minGenomes;

    /** maximum number of genomes per set */
    @Option(name = "--max", metaVar = "3000", usage = "maximum number of genomes per output set")
    private int maxGenomes;

    /** name of the input list file (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "repXXX.list.tbl", usage = "name of the input list file (if not STDIN)")
    private File inFile;


    /**
     * This utility class represents a genome from the input file.
     */
    protected class GenomeRecord implements Comparable<GenomeRecord> {

        /** genome ID */
        private final String genomeId;
        /** genome name */
        private final String genomeName;
        /** kmer score */
        private final int score;

        /**
         * Construct a genome record.
         *
         * @param line  		input line for the genome
         */
        public GenomeRecord(TabbedLineReader.Line line) {
            this.genomeId = line.get(BuildSetsProcessor.this.idCol);
            this.genomeName = line.get(BuildSetsProcessor.this.nameCol);
            this.score = line.getInt(BuildSetsProcessor.this.scoreCol);
        }

        @Override
        public int compareTo(GenomeRecord o) {
            // High score sorts first, then we sort by genome ID.
            int retVal = Integer.compare(o.score, this.score);
            if (retVal == 0)
                retVal = this.genomeId.compareTo(o.genomeId);
            return retVal;
        }

        /**
         * Write the header line for the output file.
         * 
         * @param writer    print writer for the output file
         */
        public static void writeHeader(PrintWriter writer) {
            writer.println("genome_id\tgenome_name\tscore");
        }

        /**
         * Write this genome to the output file.
         * 
         * @param writer    print writer for the output file
         */
        public void write(PrintWriter writer) {
            writer.format("%s\t%s\t%d%n", this.genomeId, this.genomeName, this.score);
        }

    }


    @Override
    protected File setDefaultOutputDir(File curDir) {
        return new File(curDir, "Sets");
    }

    @Override
    protected void setMultiReportDefaults() {
        this.minGenomes = 1500;
        this.maxGenomes = 2000;
        this.inFile = null;
    }

    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        if (this.minGenomes < 1)
            throw new ParseFailureException("Minimum number of genomes per set must be at least 1.");
        if (this.maxGenomes < this.minGenomes)
            throw new ParseFailureException("Maximum number of genomes per set must be at least the minimum.");
        // Connect to the input file.
        if (this.inFile == null) {
            log.info("Genome list will be read from standard input.");
            this.inStream = new TabbedLineReader(System.in);
        } else if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        else {
            log.info("Reading genome list from {}.", this.inFile);
            this.inStream = new TabbedLineReader(this.inFile);
        }
        // Validate the input file columns.
        this.idCol = this.inStream.findField("genome_id");
        this.nameCol = this.inStream.findField("genome_name");
        this.repCol = this.inStream.findField("rep_id");
        this.scoreCol = this.inStream.findField("score");
    }

    @Override
    protected void runMultiReports() throws Exception {
        try {
            // Create the representative map.
            this.repMap = new HashMap<>();
            // Loop through the input file, reading genomes and placing them in the map.
            log.info("Reading genome list.");
            long lastMsg = System.currentTimeMillis();
            int gCount = 0;
            for (TabbedLineReader.Line line : this.inStream) {
                GenomeRecord genome = new GenomeRecord(line);
                gCount++;
                String repId = line.get(this.repCol);
                this.repMap.computeIfAbsent(repId, k -> new ArrayList<>()).add(genome);
                // Periodically log our progress.
                long now = System.currentTimeMillis();
                if (now - lastMsg > 10000) {
                    log.info("{} genomes processed.", gCount);
                    lastMsg = now;
                }
            }
            log.info("{} genomes read for {} representatives.", gCount, this.repMap.size());
            // Now we loop through the representatives, creating output files.
            int setCount = 0;
            int skipCount = 0;
            log.info("Creating output sets.");
            for (Map.Entry<String, List<GenomeRecord>> repEntry : this.repMap.entrySet()) {
                List<GenomeRecord> genomeList = repEntry.getValue();
                if (genomeList.size() < this.minGenomes)
                    skipCount++;
                else {
                    // This representative has enough genomes. Create the output file.
                    setCount++;
                    String rep_id = repEntry.getKey();
                    String outFile = rep_id + ".tbl";
                    log.info("Writing set {} for representative {} with {} genomes.",
                            setCount, rep_id, genomeList.size());
                    try (PrintWriter writer = this.openReport(outFile)) {
                        GenomeRecord.writeHeader(writer);
                        // Sort the genome list by score.
                        Collections.sort(genomeList);
                        // Write the top N genomes.
                        int outCount = 0;
                        for (int i = 0; i < this.maxGenomes && i < genomeList.size(); i++) {
                            genomeList.get(i).write(writer);
                            outCount++;
                        }
                        log.info("Wrote set {} with {} genomes to {}.", setCount, outCount, outFile);
                    }
                }
            }
            log.info("All done. {} sets created; {} representatives skipped for lack of genomes.",
                    setCount, skipCount);
        } finally {
            this.inStream.close();
        }
    }

}
