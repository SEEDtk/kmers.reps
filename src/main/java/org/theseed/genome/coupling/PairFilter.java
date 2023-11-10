/**
 *
 */
package org.theseed.genome.coupling;

import org.theseed.basic.ParseFailureException;

/**
 * This is the base class for an object that filters potentially-coupled class pairs based on group characteristics.
 *
 * @author Bruce Parrello
 */
public abstract class PairFilter {

    // FIELDS
    /** current processor */
    protected CouplesProcessor processor;

    /**
     * Construct a new group filter.
     *
     * @param processor		parent coupling processor
     */
    public PairFilter(CouplesProcessor processor) {
        this.processor = processor;
    }

    /**
     * @return TRUE if a class pair is worth keeping.
     *
     * @param pair		class pair to check
     * @param pairData	descriptor for the genomes containing the pair
     */
    public abstract boolean isSignificant(FeatureClass.Pair pair, FeatureClass.PairData pairData);

    /**
     * @return the displayable name of this filter class
     */
    protected abstract String getName();

    /**
     * Enumeration of filter types.
     */
    public static enum Type {
        SIZE, WEIGHT, WHITELIST;

        /**
         * Create a group filter of this type.
         *
         * @param processor		parent couples processor
         *
         * @return the new group filter
         *
         * @throws ParseFailureException
         */
        public PairFilter create(CouplesProcessor processor) throws ParseFailureException {
            PairFilter retVal = null;
            switch (this) {
            case SIZE:
                retVal = new SizePairFilter(processor);
                break;
            case WEIGHT:
                retVal = new WeightPairFilter(processor);
                break;
            case WHITELIST:
                retVal = new WhitePairFilter(processor);
                break;
            }
            return retVal;
        }
    }

    @Override
    public String toString() {
        return this.getName();
    }

}
