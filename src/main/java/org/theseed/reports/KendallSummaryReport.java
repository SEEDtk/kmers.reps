/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.DoubleArray;
import org.apache.commons.math3.util.ResizableDoubleArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This error summary report computes a Kendall Tau Correlation comparing the ranking of the SSU rRNA distances and the PheS
 * distances.  The theory is that a low tau indicates an error in the SSU rRNA sequence.
 *
 * @author Bruce Parrello
 *
 */
public class KendallSummaryReport extends ErrorSummaryReport {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(KendallSummaryReport.class);
    /** map of genome IDs to genome names */
    private Map<String, String> nameMap;
    /** Kendall Tau computation engine */
    private KendallsCorrelation computer;
    /** Pearson correlation computation engine */
    private PearsonsCorrelation pComputer;
    /** map of genome IDs to correlation lists */
    private Map<String, Distances> distanceMap;
    /** output file for report */
    private File outFile;
    /** estimated number of genomes */
    private static final int ESTIMATED_GENOMES = 1000;


    public KendallSummaryReport(File summaryFile) {
        this.outFile = summaryFile;
        this.nameMap = new HashMap<String, String>(ESTIMATED_GENOMES);
        this.distanceMap = new HashMap<String, Distances>(ESTIMATED_GENOMES);
        this.computer = new KendallsCorrelation();
        this.pComputer = new PearsonsCorrelation();
    }

    @Override
    public void registerGenome(String genomeId, String genomeName) {
        this.nameMap.put(genomeId, genomeName);
    }

    @Override
    public void recordDistances(String genome1, String genome2, double protDist, double rnaDist) {
        this.recordDistancePair(genome1, protDist, rnaDist);
        this.recordDistancePair(genome2, protDist, rnaDist);
    }

    /**
     * Record a distance pair for a single genome.
     *
     * @param genomeId	ID of the genome of interest
     * @param protDist	PheS distance
     * @param rnaDist	SSU rRNA distance
     */
    public void recordDistancePair(String genomeId, double protDist, double rnaDist) {
        Distances dists = this.distanceMap.computeIfAbsent(genomeId, x -> new Distances(x));
        dists.addDistances(protDist, rnaDist);
    }

    @Override
    public void writeReport() throws IOException {
        // Sort the distances.
        SortedSet<Distances> distanceList = new TreeSet<Distances>(this.distanceMap.values());
        try (PrintWriter writer = new PrintWriter(this.outFile)) {
            writer.println("genome_id\tgenome_name\ttau_corr\tpearson");
            for (Distances dist : distanceList) {
                String genomeId = dist.getGenomeId();
                writer.format("%s\t%s\t%6.4f\t%6.4f%n", genomeId, this.nameMap.get(genomeId), dist.getCorrelation(),
                        dist.getPCorrelation());
            }
        }
    }

    /**
     * This is a utility object that holds the distance arrays for each genome.
     */
    private class Distances implements Comparable<Distances> {

        private String genomeId;
        private DoubleArray protDists;
        private DoubleArray rnaDists;
        private boolean computed;
        private double tau;
        private double pc;

        /**
         * Create a new distances object.
         *
         * @param genome	ID of the target genome
         */
        protected Distances(String genome) {
            this.protDists = createArray();
            this.rnaDists = createArray();
            this.genomeId = genome;
            this.computed = false;
        }

        /**
         * @return a resizable array for storing distances
         */
        private ResizableDoubleArray createArray() {
            return new ResizableDoubleArray(ESTIMATED_GENOMES, ESTIMATED_GENOMES, ESTIMATED_GENOMES * 1.5, ResizableDoubleArray.ExpansionMode.ADDITIVE);
        }

        /**
         * Add a pair of distances.
         *
         * @param protDist		PheS protein distance
         * @param rnaDist		SSU rRNA distance
         */
        protected void addDistances(double protDist, double rnaDist) {
            this.protDists.addElement(protDist);
            this.rnaDists.addElement(rnaDist);
        }

        /**
         * @return the kendall tau correlation
         */
        protected double getCorrelation() {
            this.checkComputed();
            return this.tau;
        }

        /**
         * @return the pearson correlation
         */
        protected double getPCorrelation() {
            this.checkComputed();
            return this.pc;
        }

        /**
         * Insure that the correlations are computed.
         */
        private void checkComputed() {
            if (! this.computed) {
                double[] prots = this.protDists.getElements();
                double[] rnas = this.rnaDists.getElements();
                this.tau = KendallSummaryReport.this.computer.correlation(prots, rnas);
                this.pc = KendallSummaryReport.this.pComputer.correlation(prots, rnas);
                this.computed = true;
            }
        }

        @Override
        public int compareTo(Distances o) {
            int retVal = Double.compare(this.getPCorrelation(), o.getPCorrelation());
            if (retVal == 0) {
                retVal = Double.compare(this.getCorrelation(), o.getCorrelation());
                if (retVal == 0)
                    retVal = this.genomeId.compareTo(o.genomeId);
            }
            return retVal;
        }

        /**
         * @return the genome ID
         */
        public String getGenomeId() {
            return this.genomeId;
        }

    }

}
