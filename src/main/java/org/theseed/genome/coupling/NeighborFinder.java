/**
 *
 */
package org.theseed.genome.coupling;

import java.util.List;

/**
 * This is the base class for the neighbor-finding algorithm.  The neighbor-finder takes as input a sorted list
 * of FeatureClass.Result objects and a starting position.  It then returns all neighbors of the specified start
 * result.
 *
 * @author Bruce Parrello
 *
 */
public abstract class NeighborFinder {

    // FIELDS
    /** parent couples processor */
    private BaseCouplingProcessor processor;
    /** maximum gap size */
    private int maxGap;

    /**
     * Construct a new neighbor finder.
     *
     * @param processor		coupling processor that will be using this finder
     */
    public NeighborFinder(BaseCouplingProcessor processor) {
        this.processor = processor;
        this.maxGap = processor.getMaxGap();
    }

    /**
     * Return the neighbors of a feature.
     *
     * @param results	sorted list of results from the current genome
     * @param pos		position in the list of the result for the feature in question
     *
     * @return a sorted list of the neighbors
     */
    public abstract List<FeatureClass.Result> getNeighbors(List<FeatureClass.Result> results, int pos);

    /**
     * Types of neighbor finders supported
     */
    public static enum Type {
        ADJACENT, CLOSE;

        /**
         * @return a neighbor finder of this type.
         *
         * @param processor		coupling processor that will be using the finder
         */
        public NeighborFinder create(BaseCouplingProcessor processor) {
            NeighborFinder retVal = null;
            switch (this) {
            case ADJACENT :
                retVal = new AdjacentNeighborFinder(processor);
                break;
            case CLOSE :
                retVal = new CloseNeighborFinder(processor);
                break;
            }
            return retVal;
        }
    }

    /**
     * @return the processor using this neighbor finder
     */
    protected BaseCouplingProcessor getProcessor() {
        return processor;
    }

    /**
     * @return the maximum allowable feature gap
     */
    protected int getMaxGap() {
        return maxGap;
    }
}
