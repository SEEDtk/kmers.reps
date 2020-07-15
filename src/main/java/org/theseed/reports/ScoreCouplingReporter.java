/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.coupling.FeatureClass;
import org.theseed.genome.coupling.FeatureClass.Pair;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.SequenceKmers;

/**
 * Produce a coupling report with scores based on group size and group distances.
 *
 * The distance is computed using genome similarity.  We locate the genome's PheS sequence and build a Kmer hash from it.  The similarity count is computed for
 * each pair of genomes in a group.  The group's similarity distance score is SUMOF(1 / (sim + 1)) over all the pairs.  The real distance score is SUMOF(jaccard-distance)
 * over all the pairs. To manage this, we have a hash from genome pairs to distances.
 *
 * @author Bruce Parrello
 *
 */
public class ScoreCouplingReporter extends CouplingReporter {

    /**
     * This is a dinky little object that contains multiple distance scores.
     */
    private static class Scores {
        /** similarity distance score */
        private double simDist;
        /** real distance score */
        private double realDist;

        /**
         * Construct a score object.
         *
         * @param sim	similarity distance
         * @param real	real distance
         */
        protected Scores(SequenceKmers first, SequenceKmers second) {
            int sim = first.similarity(second);
            this.simDist = 1.0 / (double) (sim + 1);
            this.realDist = 1.0 - (sim / (double) (first.size() + second.size() - sim));
        }

        /**
         * Construct an all-zero score object.
         */
        protected Scores() {
            this.simDist = 0.0;
            this.realDist = 0.0;
        }

        /**
         * Add another score object to this one.
         *
         * @param other		other score object to add
         */
        protected void plus(Scores other) {
            this.simDist += other.simDist;
            this.realDist += other.realDist;
        }

        /**
         * @return an output string for this score object.
         */
        public String toString() {
            return String.format("%8.4f\t%8.4f", this.simDist, this.realDist);
        }

        /**
         * @return a header string for scores
         */
        public static String header() {
            return "sim_distance\tphes_distance";
        }

    }

    // FIELDS
    /** map of genome pairs to distances */
    private Map<FeatureClass.Pair, Scores> scoreMap;
    /** map of genome IDs to protein kmer objects */
    private Map<String, ProteinKmers> genomeDb;

    /**
     * Initialize this report.
     *
     * @param output		target output stream for the report
     * @param classifier	classifier used for this report
     */
    public ScoreCouplingReporter(OutputStream output, FeatureClass classifier) {
        super(output, classifier);
        // Create the genome database and the score map.
        this.scoreMap = new HashMap<FeatureClass.Pair, Scores>(100000);
        this.genomeDb = new HashMap<String, ProteinKmers>(1000);
    }

    @Override
    public void register(Genome genome) {
        // Find the first seed protein in this genome.
        ProteinKmers prot = null;
        Iterator<Feature> iter = genome.getPegs().iterator();
        while (prot == null && iter.hasNext()) {
            Feature feat = iter.next();
            if (feat.getFunction().contentEquals("Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)")) {
                prot = new ProteinKmers(feat.getProteinTranslation());
                this.genomeDb.put(genome.getId(), prot);
            }
        }
        if (prot == null)
            throw new IllegalArgumentException("Invalid genome " + genome.toString() + " has no seed protein.");
    }

    @Override
    public void writePairLine(Pair pair, Collection<String> genomes) {
        // Get our classifier.
        FeatureClass classifier = this.getClassifier();
        // This will hold the final scores.
        Scores result = new Scores();
        // The tricky part here is we need to compute the scores for every pair of genomes in the group.
        String[] gList = new String[genomes.size()];
        gList = genomes.toArray(gList);
        for (int i = 0; i < gList.length; i++) {
            ProteinKmers protI = this.genomeDb.get(gList[i]);
            for (int j = i + 1; j < gList.length; j++) {
                ProteinKmers protJ = this.genomeDb.get(gList[j]);
                Pair gPair = classifier.new Pair(gList[i], gList[j]);
                // Find the scores for this genome pair.  If we don't know them, we compute them.
                Scores scores = this.scoreMap.computeIfAbsent(gPair, k -> new Scores(protI, protJ));
                result.plus(scores);
            }
        }
        // Write the output line.
        this.print("%s\t%d\t%s", pair.toString(), genomes.size(), result.toString());
    }

    @Override
    protected String getScoreHeadings() {
        return ("size\t" + Scores.header());
    }

}
