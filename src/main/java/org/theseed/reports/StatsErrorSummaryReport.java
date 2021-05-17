/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This error summary report displays the mean and standard deviation for the discrepancy between the SSU rRNA distance
 * and the PheS distance.  A mean significantly higher than 0.05 indicates a serious problem.
 *
 * @author Bruce Parrello
 *
 */
public class StatsErrorSummaryReport extends ErrorSummaryReport {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(StatsErrorSummaryReport.class);
    /** list of results; initially contains genome IDs and names, but eventually has the stats as well */
    private List<Result> results;
    /** map of genome IDs to error statistics */
    private Map<String, SummaryStatistics> statMap;
    /** output file */
    private File outFile;
    /** estimated number of genomes (for allocating maps) */
    private static final int ESTIMATED_GENOMES = 5000;

    /**
     * Create an error summary report output to the specified file.
     *
     * @param summaryFile	output file for report
     */
    public StatsErrorSummaryReport(File summaryFile) {
        this.outFile = summaryFile;
        this.results = new ArrayList<Result>(ESTIMATED_GENOMES);
        this.statMap = new HashMap<String, SummaryStatistics>(ESTIMATED_GENOMES);
    }

    @Override
    public void registerGenome(String genomeId, String genomeName) {
        this.results.add(new Result(genomeId, genomeName));
    }

    @Override
    public void recordDistances(String genome1, String genome2, double protDist, double rnaDist) {
        // Compute the distance error.
        double error = Math.abs(protDist - rnaDist);
        // Store it in the summary object for each genome.
        this.storeError(genome1, error);
        this.storeError(genome2, error);
    }

    /**
     * Store the specified distance error for the specified genome.
     *
     * @param genome	ID of the genome of interest
     * @param error		distance error associated with the genome
     */
    private void storeError(String genome, double error) {
        SummaryStatistics stats = this.statMap.computeIfAbsent(genome, x -> new SummaryStatistics());
        stats.addValue(error);
    }

    @Override
    public void writeReport() throws IOException {
        // Start the report.
        try (PrintWriter writer = new PrintWriter(this.outFile)) {
            log.info("Sorting results from {} genomes.", this.results.size());
            // Fill in the results.
            this.results.stream().forEach(x -> this.setResult(x));
            // Sort the filled results.  This floats the worst errors to the top.
            Collections.sort(this.results);
            log.info("Writing error summary report to {}.", this.outFile);
            // Write the header, then the data lines.
            writer.println(Result.headerLine());
            for (Result result : results)
                writer.println(result.output());
        }
    }

    /**
     * Update a result with the statistics for the report.
     *
     * @param result	result to update
     */
    private void setResult(Result result) {
        SummaryStatistics stats = this.statMap.get(result.genomeId);
        result.mean = stats.getMean();
        result.stdDev = stats.getStandardDeviation();
    }

    /**
     * This object stores the results for the report and allows us to sort them with the biggest errors to the top.
     */
    private static class Result implements Comparable<Result> {

        private String genomeId;
        private String genomeName;
        private double mean;
        private double stdDev;

        /**
         * Create a result for a genome.
         *
         * @param id		ID of the genome
         * @param name	name of the genome
         */
        public Result(String id, String name) {
            this.genomeId = id;
            this.genomeName = name;
            this.mean = 0.0;
            this.stdDev = 0.0;
        }

        @Override
        public int compareTo(Result o) {
            int retVal = Double.compare(o.mean, this.mean);
            if (retVal == 0) {
                retVal = Double.compare(o.stdDev, this.stdDev);
                if (retVal == 0)
                    retVal = this.genomeId.compareTo(o.genomeId);
            }
            return retVal;
        }

        /**
         * @return the report header
         */
        public static String headerLine() {
            return "genome_id\tgenome_name\tmean_error\tstd_dev";
        }

        /**
         * @return the report line for this result
         */
        public String output() {
            return String.format("%s\t%s\t%6.4f\t%6.4f", this.genomeId, this.genomeName, this.mean, this.stdDev);
        }

    }

}
