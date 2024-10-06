/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.coupling.CouplesProcessor;
import org.theseed.genome.coupling.FeatureClass;
import org.theseed.genome.coupling.FeatureClass.Pair;
import org.theseed.io.TabbedLineReader;


/**
 * This report displays a score indicating how often two feature classes are coupled compared to how often both
 * features occur.  Only coupled pairs from another run are examined, so it needs the output of that run to
 * determine its universe of discourse.
 *
 * @author Bruce Parrello
 *
 */
public class VerifyCouplingReporter extends CouplingReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(VerifyCouplingReporter.class);
    /** set of classifier pairs of interest */
    private Set<FeatureClass.Pair> goodPairs;
    /** number of genomes containing each pair of good classes */
    private CountMap<FeatureClass.Pair> countMap;

    /**
     * Initialize this report.
     *
     * @param output		target output stream for the report
     * @param processor		parent coupling processor for this report
     */
    public VerifyCouplingReporter(OutputStream output, CouplesProcessor processor) {
        super(output, processor);
        this.countMap = new CountMap<FeatureClass.Pair>();
        this.goodPairs = new HashSet<FeatureClass.Pair>(1000);
        // We need to fill the goodpairs set.
        FeatureClass classifier = processor.getClassifier();
        try (TabbedLineReader inStream = new TabbedLineReader(processor.getOldOutput())) {
            log.info("Reading old output.");
            for (TabbedLineReader.Line line : inStream)
                this.goodPairs.add(classifier.readPair(line));
            log.info("{} pairs selected for analysis.", this.goodPairs.size());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    protected String getScoreHeadings() {
        return "size\toccurrences\tpercent";
    }

    @Override
    public void register(Genome genome) {
        // Get the classifier.
        FeatureClass classifier = this.getClassifier();
        // Get the features of the genome.
        Collection<Feature> feats = genome.getFeatures();
        // Create a set of all the classifications found.
        Set<String> classesFound = new HashSet<String>(feats.size());
        // Loop through the features, collecting class IDs.
        for (Feature feat : feats) {
            FeatureClass.Result classes = classifier.getClasses(feat);
            if (classes != null) {
                for (String classId : classes)
                    classesFound.add(classId);
            }
        }
        // Convert the set into an array and count all the interesting pairs.
        String[] classArray = new String[classesFound.size()];
        classArray = classesFound.toArray(classArray);
        for (int i = 0; i < classArray.length; i++) {
            String classI = classArray[i];
            for (int j = i + 1; j < classArray.length; j++) {
                FeatureClass.Pair pair = classifier.new Pair(classI, classArray[j]);
                if (this.goodPairs.contains(pair))
                    this.countMap.count(pair);
            }
        }
    }

    @Override
    protected void writePairLine(Pair pair, FeatureClass.PairData genomeData) {
        // Get the number of occurrences of the pair.
        int occurs = this.countMap.getCount(pair);
        // If the pair is new, we display it with an empty percentage.
        String percent = "";
        String oString = "";
        if (this.goodPairs.contains(pair)) {
            // Here the pair is being verified, so we show its percentage.
            percent = Double.toString(genomeData.size() * 100.0 / occurs);
            // Remove it from the good-pair list so it doesn't appear in the summary.
            this.goodPairs.remove(pair);
            oString = String.format("%d", occurs);
        }
        // Format the scores.
        this.print("%s\t%d\t%s\t%s", pair.toString(), genomeData.size(), oString, percent);
    }

    @Override
    protected void summarize() {
        String percent = "0.0";
        // At the end, we spool off the pairs that weren't found but occur at least once.
        for (FeatureClass.Pair pair : this.goodPairs) {
            int occurs = this.countMap.getCount(pair);
            if (occurs > 0) {
                this.print("%s\t%d\t%d\t%s", pair.toString(), 0, occurs, percent);
            }
        }
    }

}
