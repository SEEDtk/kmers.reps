package org.theseed.sequence.seeds;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.seeds.FinderKmerBatch.ProgressTracker;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This class generates a representative set of genomes using finder-kmer technology. It 
 * performs a first pass that collects key metrics about the genome set and estimates the 
 * closeness limit to use for selecting representatives. The second pass then generates
 * representative sets based on the estimated closeness limit and one or two others as
 * determined by a fuzz factor relative to the standard deviation of the closeness limit.
 * A normal distribution curve fitting is used to generate the estimated closeness, and
 * this is reflected in the algorithm.
 * 
 * The object will read the genome IDs from a provided input stream and write the
 * representative genome list to a provided output stream. Note that a temporary
 * file is used to speed up the second pass.
 * 
 * The object is not multi-threaded. To run parallel operations, multiple instances of this 
 * class should be used. The object is serially reusable, however, meaning it can be used for
 * multiple genome sets in sequence without creating a new instance.
 * 
 */
public class FinderKmerRepGenerator {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(FinderKmerRepGenerator.class);
    /** genome ID list */
    private final Set<String> genomeIds;
    /** temporary file for storing protein sequences from the first pass */
    private final File tempFile;
    /** role map for roles to use in genome comparisons */
    private final RoleMap roleMap;
    /** ordered list of role IDs for the temporary file */
    private List<String> roleIds;
    /** connection to the BV-BRC database */
    private final P3CursorConnection p3;
    /** distance from the computed closeness limit, relative to the standard deviation, for alternate limits to try */
    private final double fuzz;
    /** batch size for processing genome IDs */
    private final int batchSize;
    /** target number of genomes for the output representative-genome set */
    private int targetSize;
    /** closeness limit for initial sample (see NOTE) */
    public static final double SAMPLE_LIMIT = 0.8;
    /** array for the above */
    private static final double[] SAMPLE_LIMIT_ARRAY = new double[] { SAMPLE_LIMIT };
    /** Z80 of the sample limit (see NOTE) */
    private static final double Z80 = 0.8416212335729143;
    /** batch size for the second pass */
    private static final int SMALL_BATCH_SIZE = 30;

    /** NOTE the Z80 value is the z-score corresponding to the 80th percentile of the standard normal distribution.
     * If the SAMPLE_LIMIT changes, the Z80 value can be computed as 
     * 
     *   import org.apache.commons.math3.distribution.NormalDistribution;
     *
     *   NormalDistribution standardNormal = new NormalDistribution();
     *   double Z80 = standardNormal.inverseCumulativeProbability(0.8); 
     */

    /**
     * Construct a new object for generating representative genomes.
     * 
     * @param roleFile      file containing the role definitions
     * @param fuzz          distance from the computed closeness limit, relative to the standard deviation, for alternate limits to try
     * @param batchSize     batch size for processing genome IDs
     * @param targetSize    target number of genomes for the output representative-genome set
     * 
     * @throws IOException  if there is an error reading the role file
     */
    public FinderKmerRepGenerator(File roleFile, double fuzz, int batchSize, int targetSize) throws IOException {
        this.roleMap = RoleMap.load(roleFile);
        log.info("{} roles found in role file {}.", this.roleMap.size(), roleFile);
        this.fuzz = fuzz;
        // Validate the fuzz factor.
        if (this.fuzz < 0.0)
            throw new IllegalArgumentException("Fuzz factor cannot be negative.");
        this.batchSize = batchSize;
        this.targetSize = targetSize;
        if (this.targetSize <= 0)
            throw new IllegalArgumentException("Target size must be positive.");
        log.info("Repgen fuzz factor is {} for target size {}.", this.fuzz, this.targetSize);
        log.info("Query batch size is {}.", this.batchSize);
        this.tempFile = File.createTempFile("repgen", ".tmp");
        log.info("Temporary file is {}.", this.tempFile);
        // Create a set to hold the genome IDs.
        this.genomeIds = new HashSet<>();
        // Connect to the BV-BRC database.
        this.p3 = new P3CursorConnection();
    }

    /**
     * Clean up the resources used by this object, such as the temporary file and the
     * genome data.
     */
    public void close() {
        if (this.tempFile != null && this.tempFile.exists()) {
            this.tempFile.delete();
        }
        if (this.genomeIds != null) {
            this.genomeIds.clear();
        }
        if (this.roleIds != null) {
            this.roleIds = null;
        }
    }

    /**
     * Set the target number of genomes for the output representative-genome set. This
     * also performs a close, resetting the internal state of the object.
     */
    public void setTargetSize(int targetSize) {
        this.close();
        this.targetSize = targetSize;
        if (this.targetSize <= 0)
            throw new IllegalArgumentException("Target size must be positive.");
    }

    /**
     * Read the genome IDs from an input streak and produce a representative genome set in the specified output stream.
     * 
     * @param inputStream   the input stream containing genome IDs
     * @param genomeIdCol   the index (1-based) or name of the genome ID column in the input stream
     * @param writer        the writer for the output representative genome set
     * 
     * @throws IOException
     */
    public void processGenomes(TabbedLineReader inputStream, String genomeIdCol, PrintWriter writer) throws IOException {
        // Read the genome IDs from the input stream.
        int colIdx = inputStream.findField(genomeIdCol);
        log.info("Reading genome IDs from column {}.", colIdx + 1);
        this.genomeIds.clear();
        inputStream.stream().forEach(line -> this.genomeIds.add(line.get(colIdx)));
        log.info("{} genome IDs found.", this.genomeIds.size());
        // Check for the edge cases.
        if (this.genomeIds.isEmpty()) {
            log.warn("No genomes to process.");
        } else if (this.genomeIds.size() <= this.targetSize) {
            log.warn("Number of genomes ({}) is less than or equal to the target size ({}).", this.genomeIds.size(), this.targetSize);
            this.outputGenomes(writer, this.genomeIds);
        } else if (this.targetSize == 1) {
            log.warn("Target size is 1, only one representative genome will be selected.");
            Iterator<String> iter = this.genomeIds.iterator();
            if (iter.hasNext())
                this.outputGenomes(writer, Collections.singleton(iter.next()));
        }
        log.info("Fuzz factor is {}", this.fuzz);
        // For the first pass, we will loop through the genome IDs, building and processing
        // batches. Each batch will be used to update a statistics object for the closeness values
        // and create a repgen set for the sample closeness limit. The size of that repgen set
        // will be one of the key metrics for determining the final closeness limit. Finally, we
        // save the ID, name, and protein sequences for each genome in a temporary file so we can
        // make the second pass without using the database.
        log.info("Preparing for first pass through {} genomes.", this.genomeIds.size());
        // Create a progress tracker so we can monitor the processing of genome batches.
        ProgressTracker progress = new ProgressTracker();
        // Create a statistics object to track the closeness values and attach it to the progress tracker.
        FinderKmerStats stats = new FinderKmerStats();
        progress.setStats(stats);
        // This will be used to hold the size of the sample repgen set.
        int sampleSize;
        // Open the temporary file for output.
        try (PrintWriter tempWriter = new PrintWriter(this.tempFile)) {
            // Create a FinderKmerBatch to hold the sample repgen set.
            FinderKmerBatch repgen = new FinderKmerBatch();
            // Put it into an array because of the way the processBatch method expects it.
            FinderKmerBatch[] repgenArray = new FinderKmerBatch[] { repgen };
            // Write the header for the temporary file.
            this.roleIds = FinderKmerBatch.createHeader(roleMap, tempWriter);
            // We will loop through the genome IDs, filling the sample batch. Once its full,
            // we generate the repgen set, update the statistics, and cache the proteins in
            // the temporary file.
            Set<String> batch = new HashSet<>();
            for (String genomeId : this.genomeIds) {
                batch.add(genomeId);
                if (batch.size() >= this.batchSize) {
                    this.processBatch(batch, tempWriter, repgenArray, progress);
                    batch.clear();
                }
            }
            // Process any remaining genomes in the batch.
            if (! batch.isEmpty()) {
                this.processBatch(batch, tempWriter, repgenArray, progress);
                batch.clear();
            }
            // Record the size of the sample repgen set.
            sampleSize = repgen.size();
        }
        double desiredCloseness = estimateY(stats, this.genomeIds.size(), sampleSize, (double) this.targetSize);
        log.info("Estimated desired closeness: {}", desiredCloseness);
        // We are going to make three repgen sets with three closeness limits. First, we build the array of closeness limits.
        double fuzzValue = stats.getStdDev() * this.fuzz;
        log.info("Fuzz value is {}", fuzzValue);
        double[] closenessLimits = this.computeClosenessLimits(desiredCloseness, fuzzValue);
        FinderKmerBatch[] repgenArray = this.buildRepgenArray(closenessLimits);
        // Reset the progress tracker for the second pass.
        progress.reset();
        log.info("Computing repgen sets for the second pass.");
        // Prepare to read back the temporary file.
        try (TabbedLineReader tempReader = new TabbedLineReader(this.tempFile)) {
            this.roleIds = FinderKmerBatch.loadHeader(tempReader);
            // This object will hold the current batch of genomes.
            FinderKmerBatch batch = new FinderKmerBatch();
            // Count the number of lines read.
            int lineCount = 0;
            while (tempReader.hasNext()) {
                batch.load(roleIds, SMALL_BATCH_SIZE, tempReader);
                lineCount += batch.size();
                batch.updateRepgens(repgenArray, closenessLimits, progress);
                batch.clear();
            }
            log.info("Completed second pass with {} lines processed.", lineCount);
        }
        FinderKmerBatch repgen = this.selectBestRepgen(repgenArray);
        log.info("Selected best repgen set with {} genomes.", repgen.size());
        // Output the representative genomes.
        this.outputGenomes(writer, repgen.genomeIds());
    }

    /**
     * Process a batch of genomes during the first pass. We load the genomes from the database, update the statistics
     * merge them into the sample repgen set, and cache the proteins in the temporary file.
     *  
     * @param batch         the set of genome IDs to process in this batch
     * @param tempWriter    the writer for the temporary file to cache proteins
     * @param repgenSample  the sample representative genome set to merge into
     * @param progress      the progress tracker for updating statistics
     * 
     * @throws IOException 
     */
    private void processBatch(Set<String> batch, PrintWriter tempWriter, FinderKmerBatch[] repgenSample,
            ProgressTracker progress) throws IOException {
        log.info("Loading batch of {} genomes for first pass.", batch.size());
        FinderKmerBatch kmerBatch = new FinderKmerBatch();
        kmerBatch.loadBatch(roleMap, p3, batch);
        // Write the cached proteins to the temporary file.
        kmerBatch.write(this.roleIds, tempWriter);
        // Merge the batch into the sample representative genome set.
        kmerBatch.updateRepgens(repgenSample, SAMPLE_LIMIT_ARRAY, progress);
    }


    /**
     * Compute the estimated closeness limit (y) for a sample with the specified statistical characteristics, given
     * a desired output repgen size (x)
     * 
     * @param stats         the statistical characteristics of the sample
     * @param size          the number of genomes in the sample
     * @param sample        repgen size at closeness 0.8
     * @param x             desired output repgen size
     * 
     * @return              the estimated closeness limit (y) for the given sample and desired output repgen size
     */
    private static double estimateY(FinderKmerStats stats, int size, double sample, double x) {
        if (sample <= 0 || size <= 0) {
            throw new IllegalArgumentException("size and sample must be positive");
        }
        if (x == sample) {
            return 0.8;
        }

        double logSigma = -1.7128
                - 1.0084 * stats.getMean()
                + 1.6362 * stats.getStdDev()
                + 0.4290 * Math.log((double) size)
                - 0.3801 * Math.log(sample);
        double sigma = Math.exp(logSigma);

        double ratio = Math.max(x, 0.5) / sample;
        double z = Math.log(ratio) / sigma + Z80;
        return normalCdf(z);
    }

    /** Standard normal CDF, Phi(z). */
    private static double normalCdf(double z) {
        return 0.5 * (1.0 + erf(z / Math.sqrt(2.0)));
    }

    /** Abramowitz and Stegun approximation 7.1.26, max error ~1.5e-7. */
    private static double erf(double x) {
        double sign = Math.signum(x);
        double ax = Math.abs(x);
        double t = 1.0 / (1.0 + 0.3275911 * ax);
        double y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
                - 0.284496736) * t + 0.254829592) * t * Math.exp(-ax * ax);
        return sign * y;
    }    

    /**
     * Compute the closeness limits to use for selecting representative genomes. The base value is
     * the desired closeness. We will also attempt to use the desired closeness plus or minus the fuzz value.
     * If a computed closeness is outside the acceptable range (0.0 to 1.0), it will be ignored.
     * 
     * @param desiredCloseness    the base closeness value to use
     * @param fuzzValue           the amount to vary the closeness above and below the base value
     * @return                    an array of closeness limits to try, including the base value and the fuzzed values
     */
    private double[] computeClosenessLimits(double desiredCloseness, double fuzzValue) {
        List<Double> limits = new ArrayList<>();
        double base = desiredCloseness;
        if (base >= 0.0 && base <= 1.0) {
            limits.add(base);
        }
        // Insure the user has not turned off the fuzz factor.
        if (fuzzValue > 0.0) {
            double lower = base - fuzzValue;
            double upper = base + fuzzValue;
            if (lower >= 0.0) {
                limits.add(lower);
            }
            if (upper <= 1.0) {
                limits.add(upper);
            }
        }
        log.info("Computed closeness limits: {}", StringUtils.join(limits, ", "));
        return limits.stream().mapToDouble(Double::doubleValue).toArray();
    }

    /**
     * Build an array of empty FinderKmerBatch objects, one per closeness limit provided,
     * to hold the repgen sets.
     * 
     * @param closenessLimits    the array of closeness limits for which to create empty FinderKmerBatch objects
     * @return                    an array of empty FinderKmerBatch objects, one per closeness limit
     */
    private FinderKmerBatch[] buildRepgenArray(double[] closenessLimits) {
        FinderKmerBatch[] retVal = new FinderKmerBatch[closenessLimits.length];
        for (int i = 0; i < closenessLimits.length; i++) {
            retVal[i] = new FinderKmerBatch();
        }
        return retVal;
    }

    /**
     * Find the repgen set whose size is closest to the desired target size.
     * 
     * @param repgenArray    the array of FinderKmerBatch objects to evaluate
     * 
     * @return               the FinderKmerBatch object whose size is closest to the desired target size
     */
    private FinderKmerBatch selectBestRepgen(FinderKmerBatch[] repgenArray) {
        FinderKmerBatch retVal = null;
        int error = Integer.MAX_VALUE;
        for (FinderKmerBatch batch : repgenArray) {
            int currentError = Math.abs(batch.size() - this.targetSize);
            if (currentError < error) {
                error = currentError;
                retVal = batch;
            }
        }
        return retVal;
    }

    /**
     * Write the representative genomes to the output writer. This involves looking up the genome names
     * from the database.
     * 
     * @param writer        the output stream for the representative genome information
     * @param repIds        the collection of representative genome IDs to output
     * 
     * @throws IOException 
     */
    private void outputGenomes(PrintWriter writer, Collection<String> repIds) throws IOException {
        log.info("Retrieving genome names.");
        // Get the genome names.
        List<JsonObject> genomes = this.p3.getRecords("genome", P3CursorConnection.MAX_LIMIT, this.batchSize, "genome_id",
                repIds, "genome_id,genome_name");
        log.info("Writing {} genomes to output.", genomes.size());
        writer.println("genome_id\tgenome_name");
        for (JsonObject genome : genomes) {
            String genomeId = KeyBuffer.getString(genome, "genome_id");
            String genomeName = KeyBuffer.getString(genome, "genome_name");
            writer.println(genomeId + "\t" + genomeName);
        }
        writer.flush();
    }

}
