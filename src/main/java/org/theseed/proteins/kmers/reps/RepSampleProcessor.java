package org.theseed.proteins.kmers.reps;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.kohsuke.args4j.Option;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.sequence.seeds.FinderKmerBatch;
import org.theseed.utils.BasePipeProcessor;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This command builds a representative set of genomes using a finderKmer load file. This provides a balanced sample
 * of the genome space for study and analysis. The selection is based on a similarity (closeness) threshold to ensure 
 * diversity among the representatives. The finderkmer load file contains protein sequences used to determine the
 * closeness.
 * 
 * The finderKmer load file will be read from the standard input. The standard output will be the IDs and names
 * of a representative set of genomes. The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -i   input finderKmer load file (if not STDIN)
 * -o   output file for the representative genome set (if not STDOUT)
 * -b   batch size for database queries (default is 400)
 * 
 * --sim    similarity threshold for representative selection (default is 0.9)
 * 
 * RepSampleProcessor
 */
public class RepSampleProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(RepSampleProcessor.class);
    /** finder-kmer batch containing current representative genomes */
    private FinderKmerBatch repBatch;
    /** finder-kmer batch containing new genomes to check */
    private FinderKmerBatch newBatch;
    /** list of role IDs from the load file */
    private List<String> roleIds;
    /** batch size for loading genomes from the load file */
    private static final int LOAD_BATCH_SIZE = 30;

    // COMMAND-LINE OPTIONS

    /** batch size for database queries */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "100", usage = "batch size for database queries")
    private int batchSize;

    /** similarity threshold for representative selection */
    @Option(name = "--sim", metaVar = "0.95", usage = "similarity threshold for representative selection")
    private double simThreshold;


    @Override
    protected void setPipeDefaults() {
        this.batchSize = 400;
        this.simThreshold = 0.9;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Verify that the first column is genome IDs, and then get the list of role IDs from the header.
        String[] headers = inputStream.getLabels();
        if (! headers[0].equals("genome_id"))
            throw new IOException("First column must be genome IDs");
        this.roleIds = Arrays.asList(Arrays.copyOfRange(headers, 1, headers.length));
        log.info("{} role IDs found in the load file.", this.roleIds.size());
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Insure we have a valid batch size and similarity threshold.
        if (this.batchSize <= 0)
            throw new IOException("Batch size must be positive");
        if (this.simThreshold <= 0.0 || this.simThreshold > 1.0)
            throw new IOException("Similarity threshold must be between 0 and 1");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Create the finderkmer batch for the representative genomes.
        this.repBatch = new FinderKmerBatch();
        // Create the finderkmer batch for the new genomes from the load file.
        this.newBatch = new FinderKmerBatch();
        // Count the genomes read and batches processed.
        int genomeCount = 0;
        int batchCount = 0;
        long startTime = System.currentTimeMillis();
        while (inputStream.hasNext()) {
            this.newBatch.load(this.roleIds, LOAD_BATCH_SIZE, inputStream);
            genomeCount += this.newBatch.size();
            batchCount++;
            log.info("Processing batch {} with {} genomes.", batchCount, this.newBatch.size());
            // Process the current batch and clear it.
            this.repBatch.processBatch(this.newBatch, this.simThreshold);
            this.newBatch.clear();
            long elapsedTime = System.currentTimeMillis() - startTime;
            int ngPerMin = (int) (60000.0 * genomeCount / elapsedTime);
            log.info("{} representatives found out of {} genomes ({} per minute).", this.repBatch.size(), genomeCount,
                    ngPerMin);
            
        }
        log.info("Processing complete. {} representatives found out of {} genomes.", this.repBatch.size(), genomeCount);
        // Now we must write out the genomes. First, we need the genome names.
        // Get the genome names for the representatives.
        Map<String,String> genomeMap = this.getGenomeNames();
        log.info("{} genome names found.", genomeMap.size());
        // Write the output header.
        writer.println("genome_id\tgenome_name");
        // Write the representative genomes.
        for (String genomeId : this.repBatch.genomeIds()) {
            String genomeName = genomeMap.getOrDefault(genomeId, "");
            writer.println(genomeId + "\t" + genomeName);
        }
        log.info("{} representative genomes written to output.", this.repBatch.size());
    }

    /**
     * Read the names of the representative genomes from the BV-BRC database.
     * 
     * @return a map of genome IDs to genome names for the representative genomes
     * 
     * @throws IOException if there is an error accessing the database
     */
    private Map<String, String> getGenomeNames() throws IOException {
        // Connect to the BV-BRC database.
        P3CursorConnection p3 = new P3CursorConnection();
        // Get the genome IDs from the representative-genome batch.
        Set<String> repGenomeIds = this.repBatch.genomeIds();
        // Query the database for the genome names.
        List<JsonObject> repGenomes = p3.getRecords("genome", P3CursorConnection.MAX_LIMIT, this.batchSize,
                "genome_id", repGenomeIds, "genome_id,genome_name");
        // Convert the list of JSON objects to a map of genome IDs to genome names.
        Map<String, String> retVal = new HashMap<>(repGenomes.size() * 3 / 2 + 1);
        for (JsonObject genome : repGenomes) {
            String genomeId = KeyBuffer.getString(genome, "genome_id");
            String genomeName = KeyBuffer.getString(genome, "genome_name");
            retVal.put(genomeId, genomeName);
        }
        return retVal;
    }

}
