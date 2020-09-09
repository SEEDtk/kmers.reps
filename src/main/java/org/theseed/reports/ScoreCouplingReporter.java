/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.coupling.BaseCouplingProcessor;
import org.theseed.genome.coupling.FeatureClass;
import org.theseed.proteins.Function;
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
        /** class1 percentage */
        private double percent1;
        /** class2 percentage */
        private double percent2;
        /** class1 size */
        private int size1;
        /** class2 size */
        private int size2;

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
         * Finish the scores.
         *
         * @param pair		pair being scored
         * @param size		size of the group
         */
        protected void finish(ScoreCouplingReporter parent, FeatureClass.Pair pair, int size) {
            double numerator = size * 100;
            this.size1 = parent.classCounts.getCount(pair.getClass1());
            this.size2 = parent.classCounts.getCount(pair.getClass2());
            this.percent1 = numerator / this.size1;
            this.percent2 = numerator / this.size2;
        }

        /**
         * @return an output string for this score object.
         */
        public String toString() {
            return String.format("%8.4f\t%8.4f\t%8.2f\t%8.2f\t%d\t%d", this.simDist, this.realDist,
                    this.percent1, this.percent2, this.size1, this.size2);
        }

        /**
         * @return a header string for scores
         */
        public static String header() {
            return "sim_distance\tphes_distance\tpercent1\tpercent2\tfamily1\tfamily2";
        }

    }

    // FIELDS
    /** map of genome pairs to distances */
    private Map<FeatureClass.Pair, Scores> scoreMap;
    /** map of genome IDs to protein kmer objects */
    private Map<String, ProteinKmers> genomeDb;
    /** count of genomes containing each class */
    private CountMap<String> classCounts;
    /** functional assignment for the seed protein */
    private static final String SEED_PROTEIN_FUNCTION = "Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)";


    /**
     * Initialize this report.
     *
     * @param output		target output stream for the report
     * @param processor		parent coupling processor for this report
     */
    public ScoreCouplingReporter(OutputStream output, BaseCouplingProcessor processor) {
        super(output, processor);
        // Create the genome database and the score map.
        this.scoreMap = new HashMap<FeatureClass.Pair, Scores>(100000);
        this.genomeDb = new HashMap<String, ProteinKmers>(1000);
        this.classCounts = new CountMap<String>();
    }

    @Override
    public void register(Genome genome) {
        // We need to count the classes that occur in this genome and find the seed protein.
        FeatureClass classifier = this.getClassifier();
        // First, we get the features.
        Collection<Feature> feats = genome.getFeatures();
        // This will contain the longest seed protein.
        String prot = "";
        // This will track the features found.
        Set<String> classes = new HashSet<String>();
        // Loop through the features.
        for (Feature feat : feats) {
            // Check for the seed protein.  We keep the longest.
            String function = feat.getFunction();
            if (function != null && Function.commentFree(function).contentEquals(SEED_PROTEIN_FUNCTION)) {
                String newProt = feat.getProteinTranslation();
                if (newProt != null && newProt.length() > prot.length()) prot = newProt;
            }
            // Check for the classes.
            FeatureClass.Result result = classifier.getClasses(feat);
            if (result != null) {
                for (String classId : result)
                    classes.add(classId);
            }
        }
        // Store the seed protein.
        if (prot.isEmpty())
            throw new IllegalArgumentException("Invalid genome " + genome.toString() + " has no seed protein.");
        else
            this.genomeDb.put(genome.getId(), new ProteinKmers(prot));
        // Count the classes.
        for (String classId : classes)
            this.classCounts.count(classId);
    }

    @Override
    public void writePairLine(FeatureClass.Pair pair, FeatureClass.PairData genomeData) {
        // Get our classifier.
        FeatureClass classifier = this.getClassifier();
        // This will hold the final scores.
        Scores result = new Scores();
        // The tricky part here is we need to compute the scores for every pair of genomes in the group.
        String[] gList = new String[genomeData.size()];
        gList = genomeData.getGenomes().toArray(gList);
        for (int i = 0; i < gList.length; i++) {
            ProteinKmers protI = this.genomeDb.get(gList[i]);
            for (int j = i + 1; j < gList.length; j++) {
                ProteinKmers protJ = this.genomeDb.get(gList[j]);
                FeatureClass.Pair gPair = classifier.new Pair(gList[i], gList[j]);
                // Find the scores for this genome pair.  If we don't know them, we compute them.
                Scores scores = this.scoreMap.computeIfAbsent(gPair, k -> new Scores(protI, protJ));
                result.plus(scores);
            }
        }
        // Finish off the scores.
        result.finish(this, pair, genomeData.size());
        // Write the output line.
        this.print("%s\t%d\t%8.4f\t%s", pair.toString(), genomeData.size(), genomeData.weight(), result.toString());
    }

    @Override
    protected String getScoreHeadings() {
        return ("size\tweight\t" + Scores.header());
    }

}
