/**
 *
 */
package org.theseed.genome.coupling;

import org.theseed.genome.coupling.FeatureClass.Pair;
import org.theseed.genome.coupling.FeatureClass.PairData;
import org.theseed.utils.ParseFailureException;

/**
 * This pair filter only accepts pairings whose genome group has a weighted size at or above a certain limit.
 * The weight of a genome is inversely proportional to the product of how many times each class occurs in
 * the genome.
 *
 * @author Bruce Parrello
 *
 */
public class WeightPairFilter extends PairFilter {

    // FIELDS
    /** minimum group weight */
    private double minWeight;

    public WeightPairFilter(CouplesProcessor processor) throws ParseFailureException {
        super(processor);
        this.minWeight = processor.getGroupLimit();
        if (this.minWeight < 0.0)
            throw new ParseFailureException("Minimum weight limit must be non-negative.");
    }

    @Override
    public boolean isSignificant(Pair pair, PairData pairData) {
        return (pairData.weight() >= this.minWeight);
    }

    @Override
    protected String getName() {
        return "WEIGHT";
    }

}
