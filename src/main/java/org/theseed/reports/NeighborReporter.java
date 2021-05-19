/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.CloseSets;

/**
 * This is the base class for reports produced by the BadRnaProcessor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class NeighborReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(NeighborReporter.class);
    /** output writer */
    private PrintWriter writer;

    /**
     * This enum describes the report types.
     */
    public static enum Type {
        STATS {
            @Override
            public NeighborReporter create(PrintWriter writer) {
                return new StatsNeighborReporter(writer);
            }
        }, SETS {
            @Override
            public NeighborReporter create(PrintWriter writer) {
                return new SetsNeighborReporter(writer);
            }
        };

        /**
         * @return a neighborhood reporter of this type
         *
         * @param writer	print writer to receive the report output
         */
        public abstract NeighborReporter create(PrintWriter writer);
    }
    /**
     * Construct a new neighbor report.
     *
     * @param writer	output writer for report
     */
    public NeighborReporter(PrintWriter writer) {
        this.writer = writer;
    }

    /**
     * Write a formatted output line.
     */
    protected void print(String format, Object... args) {
        this.writer.format(format, args);
        this.writer.println();
    }

    /**
     * Write an unformatted output line.
     */
    protected void println(String line) {
        this.writer.println(line);
    }

    /**
     * Write a blank output line.
     */
    protected void println() {
        this.writer.println();
    }

    /**
     * Initialize the report.
     */
    public abstract void openReport();

    /**
     * Record the testing results for a genome.
     *
     * @param id		ID of the genome
     * @param name		name of the genome
     * @param results	distance results for the genome against the reference set
     */
    public abstract void recordGenome(String id, String name, CloseSets results);

    /**
     * Close and finish the report.
     */
    public abstract void finish();

}
