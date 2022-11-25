/**
 *
 */
package org.theseed.genome.coupling;

import java.util.Collection;

import org.theseed.counters.Shuffler;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader.Line;

/**
 * This is a special classification method used to test the efficacy of various parameters.  It
 * assigns a randomly-selected feature class ID to each peg.
 *
 * @author Bruce Parrello
 *
 */
public class RandomFeatureClass extends FeatureClass {

    // FIELDS
    /** current genome's ID */
    private String genomeId;
    /** set of protein families in this genome */
    private Shuffler<String> families;
    /** current position in the shuffler */
    private int pos;

    /**
     * Construct this object.
     */
    public RandomFeatureClass() {
        this.genomeId = "";
        this.families = new Shuffler<String>(1000);
        this.pos = 0;
    }

    @Override
    public Result getClasses(Feature feat) {
        Result retVal = new Result(feat);
        if (feat.getPgfam() != null) {
            this.validateGenome(feat.getParent());
            String classId = families.get(this.pos);
            this.pos++;
            retVal.add(classId);
        }
        return retVal;
    }

    /**
     * Insure we have the list of families in the current genome available.
     *
     * @param genome	genome for the current feature.
     */
    private void validateGenome(Genome genome) {
        if (! this.genomeId.contentEquals(genome.getId())) {
            // Here we need to refresh the list.
            this.families.clear();
            for (Feature peg : genome.getPegs()) {
                String family = peg.getPgfam();
                if (family != null)
                    families.add(family);
            }
            // Shuffle the families.
            this.families.shuffle(families.size());
            // Position on the first.
            this.pos = 0;
        }
    }

    @Override
    public String getName(String classId) {
        return classId;
    }

    @Override
    public String getHeadings() {
        return "class1\tclass2";
    }

    @Override
    public void cacheNames(Collection<String> classes) { }

    @Override
    public Pair readPair(Line line) {
        return new Pair(line.get(0), line.get(1));
    }

}
