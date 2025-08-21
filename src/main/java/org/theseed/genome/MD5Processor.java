/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.LinkedHashMap;
import java.util.Map;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3MD5Hex;

/**
 * This command reads a tab-delimited file containing genome IDs and adds the DNA MD5 checksum at
 * the end.  The tab-delimited file should be on the standard input.  The output file is written
 * to the standard output, but progress messages are sent to the log, so the log should be configured
 * to the standard error output or a file.
 *
 * The following command-line options are supported.
 *
 * -h	display command help
 * -v	log more detailed progress messages
 * -c	index (1-based) or name of the column containing genome IDs
 * -b	number of genomes to process in each batch
 * -i	the input file (default is STDIN)
 * -r	resume after an error
 *
 * @author Bruce Parrello
 *
 */
public class MD5Processor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(MD5Processor.class);
    /** input column index */
    private int gCol;
    /** input file stream */
    private TabbedLineReader inStream;
    /** current batch of input lines, keyed by genome ID */
    private LinkedHashMap<String, String> batch;
    /** MD5 computer */
    private P3MD5Hex mdComputer;
    /** ID of the last genome processed */
    private String lastGenome;

    // COMMAND-LINE OPTIONS

    /** genome ID column index */
    @Option(name = "-c", aliases = { "--col" }, metaVar = "genome_id", usage = "index (1-based) or name of input column")
    private String gColumn;

    /** number of genomes to process at a time */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "20", usage = "number of genomes per batch")
    private int batchSize;

    @Option(name = "-i", aliases = { "--input" }, metaVar = "inFile", usage = "input file (if not STDIN)")
    private File inFile;

    @Option(name = "-r", aliases = { "--resume" }, metaVar = "562.36375", usage = "resume after the specified genome")
    private String resume;

    @Override
    protected void setDefaults() {
        this.gColumn = "genome_id";
        this.batchSize = 10;
        this.inFile = null;
        this.resume = null;
        this.lastGenome = "none";
    }

    @Override
    protected void validateParms() throws IOException {
        // Set up the input stream.
        if (this.inFile == null) {
            this.inStream = new TabbedLineReader(System.in);
        } else if (! this.inFile.canRead()) {
            throw new FileNotFoundException("Input file " + this.inFile + " not found or unreadable.");
        } else {
            this.inStream = new TabbedLineReader(this.inFile);
        }
        // Find the input column.
        this.gCol = this.inStream.findField(this.gColumn);
    }

    @Override
    public void runCommand() {
        try {
            log.info("Preparing for processing.");
            this.mdComputer = new P3MD5Hex();
            // The basic approach is to read in a batch of genome IDs.  We accumulate input lines until we have a batch
            // and then compute the MD5s.
            this.batch = new LinkedHashMap<>(this.batchSize);
            // If we are not resuming, we write the output header.
            if (this.resume == null)
                System.out.println(this.inStream.header() + "\tmd5_checksum");
            long genomeCount = 0;
            // For resume mode, skip to the specified genome.
            if (this.resume != null) {
                log.info("Skipping ahead past {}.", this.resume);
                TabbedLineReader.Line line = this.inStream.next();
                genomeCount++;
                while (! line.get(this.gCol).contentEquals(this.resume)) {
                    line = this.inStream.next();
                    genomeCount++;
                }
                log.info("{} lines skipped.", genomeCount);
                this.lastGenome = this.resume;
            }
            // Now loop through the data lines.
            long start = System.currentTimeMillis();
            long processCount = 0;
            for (TabbedLineReader.Line line : this.inStream) {
                // Insure there is room for another genome in this batch.
                if (this.batch.size() >= this.batchSize) {
                    this.processBatch();
                    this.batch.clear();
                    log.info("{} genomes processed. {} seconds/genome.", genomeCount,
                            ((double) (System.currentTimeMillis() - start) / (1000 * processCount)));
                }
                // Add this genome to the batch.
                String genomeId = line.get(this.gCol);
                this.batch.put(genomeId, line.getAll());
                genomeCount++;
                processCount++;
            }
            // Process the residual batch.
            this.processBatch();
            log.info("All done. {} genomes processed.", genomeCount);
        } catch (UnsupportedEncodingException | NoSuchAlgorithmException e) {
            log.error("Error after genome {}.", this.lastGenome, e);
        } finally {
            this.inStream.close();
        }
    }

    /**
     * Process a batch of genomes and output the data lines with their MD5s.
     *
     * @throws UnsupportedEncodingException
     */
    private void processBatch() throws UnsupportedEncodingException {
        Map<String,String> md5Map = this.mdComputer.genomeMD5s(this.batch.keySet());
        // Now we loop through the genomes in the order they were read (a characteristic of the LinkedHashMap.
        for (String genome : this.batch.keySet()) {
            String md5 = md5Map.get(genome);
            if (md5 == null)
                log.info("{} not found in database or has no contigs.", genome);
            else {
                System.out.println(this.batch.get(genome) + "\t" + md5);
                this.lastGenome = genome;
            }
        }
    }

}
