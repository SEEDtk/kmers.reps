package org.theseed.proteins.kmers.reps;

/**
 * This interface is used for objects that store a list of representative-genome databases.
 *
 * @author Bruce Parrello
 *
 */
public interface IRepGenContainer {

    /**
     * Create the array of RepGen databases.
     *
     * @param size	estimated number needed
     */
    void initRepGenSets(int size);

    /**
     * Store a new RepGen database in the list.
     *
     * @param repGenomeDb	RepGen database to store
     */
    void addRepGenSet(RepGenomeDb repGenomeDb);

}