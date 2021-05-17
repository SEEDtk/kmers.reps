/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;

/**
 * This object manages the various types of error summary report.  Currently, only two are supported:  a null
 * report that does nothing, and a full report that is output to a file.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ErrorSummaryReport {

    /**
     * Save the name of a genome.
     *
     * @param genomeId		ID of the genome
     * @param genomeName	name of the genome
     */
    public abstract void registerGenome(String genomeId, String genomeName);

    /**
     * Record genome distances.
     *
     * @param genome1	ID of the first genome
     * @param genome2	ID of the second genome
     * @param protDist	seed protein distance
     * @param rnaDist	SSU rRNA distance
     */
    public abstract void recordDistances(String genome1, String genome2, double protDist, double rnaDist);

    /**
     * Produce the report.
     *
     * @throws IOException
     */
    public abstract void writeReport() throws IOException;

    /**
     * The null report saves no data and produces no output.
     */
    public static class Null extends ErrorSummaryReport {

        @Override
        public void registerGenome(String genomeId, String genomeName) {
        }

        @Override
        public void recordDistances(String genome1, String genome2, double protDist, double rnaDist) {
        }

        @Override
        public void writeReport() {
        }

    }

    /**
     * Enumeration for the type of report.
     */
    public static enum Type {
        STATS {
            @Override
            public ErrorSummaryReport create(File summaryFile) {
                return new StatsErrorSummaryReport(summaryFile);
            }
        }, KENDALL {
            @Override
            public ErrorSummaryReport create(File summaryFile) {
                return new KendallSummaryReport(summaryFile);
            }
        };

        public abstract ErrorSummaryReport create(File summaryFile);
    }

    /**
     * @return a summary reporter for the type of file specified
     *
     * @param summaryFile	the output file for the report (if NULL, a null reporter will be returned
     */
    public static ErrorSummaryReport create(Type reportType, File summaryFile) {
        ErrorSummaryReport retVal;
        if (summaryFile == null)
            retVal = new ErrorSummaryReport.Null();
        else
            retVal = reportType.create(summaryFile);
        return retVal;
    }
}
