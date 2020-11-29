/**
 *
 */
package org.theseed.genome.coupling;

import org.theseed.genome.coupling.FeatureClass.Pair;
import org.theseed.genome.coupling.FeatureClass.PairData;
import org.theseed.utils.ParseFailureException;

/**
 * @author Bruce Parrello
 *
 */
public class SizePairFilter extends PairFilter {

    // FIELDS
    /** minimum group size */
    private int minSize;

    public SizePairFilter(CouplesProcessor processor) throws ParseFailureException {
        super(processor);
        // The idea here is that 15.0 becomes 15, but 15.1 becomes 16.
        double limit = Math.ceil(processor.getGroupLimit());
        this.minSize = (int) limit;
        if (limit <= 1.0)
            throw new ParseFailureException("Minimum group size limit must be greater than 1.");
    }

    @Override
    public boolean isSignificant(Pair pair, PairData pairData) {
        return (pairData.size() >= this.minSize);
    }

    @Override
    protected String getName() {
        return "SIZE";
    }

}
