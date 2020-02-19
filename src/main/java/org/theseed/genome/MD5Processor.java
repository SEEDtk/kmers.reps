/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.LinkedHashMap;
import java.util.Map;

import org.kohsuke.args4j.Option;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3MD5Hex;
import org.theseed.utils.BaseProcessor;

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
 *
 * @author Bruce Parrello
 *
 */
public class MD5Processor extends BaseProcessor {

    // FIELDS
    /** input column index */
    private int gCol;
    /** input file stream */
    private TabbedLineReader inStream;
    /** current batch of input lines, keyed by genome ID */
    private LinkedHashMap<String, String> batch;
    /** MD5 computer */
    private P3MD5Hex mdComputer;

    // COMMAND-LINE OPTIONS

    /** genome ID column index */
    @Option(name = "-c", aliases = { "--col" }, metaVar = "genome_id", usage = "index (1-based) or name of input column")
    private String gColumn;

    /** number of genomes to process at a time */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "20", usage = "number of genomes per batch")
    private int batchSize;

    @Option(name = "-i", aliases = { "--input" }, metaVar = "inFile", usage = "input file (if not STDIN)")
    private File inFile;

    @Override
    protected void setDefaults() {
        this.gColumn = "genome_id";
        this.batchSize = 10;
        this.inFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
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
        return true;
    }

    @Override
    public void run() {
        try {
            this.mdComputer = new P3MD5Hex();
            // The basic approach is to read in a batch of genome IDs.  We accumulate input lines until we have a batch
            // and then compute the MD5s.
            this.batch = new LinkedHashMap<String, String>(this.batchSize);
            // First, we write the output header.
            System.out.println(this.inStream.header() + "\tmd5_checksum");
            // Now loop through the data lines.
            long genomeCount = 0;
            long start = System.currentTimeMillis();
            for (TabbedLineReader.Line line : this.inStream) {
                // Insure there is room for another genome in this batch.
                if (this.batch.size() >= this.batchSize) {
                    this.processBatch();
                    this.batch.clear();
                    log.info("{} genomes processed. {} seconds/genome.", genomeCount,
                            ((double) (System.currentTimeMillis() - start) / (1000 * genomeCount)));
                }
                // Add this genome to the batch.
                String genomeId = line.get(this.gCol);
                this.batch.put(genomeId, line.getAll());
                genomeCount++;
            }
            // Process the residual batch.
            this.processBatch();
            log.info("All done. {} genomes processed.", genomeCount);
        } catch (Exception e) {
            e.printStackTrace();
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
            else
                System.out.println(this.batch.get(genome) + "\t" + md5);
        }
    }

}
