/**
 *
 */
package org.theseed.genome;

import java.util.ArrayList;
import java.util.List;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * This object contains two lists of close genomes for a given target genome.
 *
 * @author Bruce Parrello
 */
public class CloseSets {

    // FIELDS
    /** set of neighbors computed using PheS */
    private NavigableSet<CloseGenome> protNeighbors;
    /** set of neighbors computed using SSU rRNA */
    private NavigableSet<CloseGenome> rnaNeighbors;

    /**
     * Create an empty neighbor set.
     */
    public CloseSets() {
        this.protNeighbors = new TreeSet<CloseGenome>();
        this.rnaNeighbors = new TreeSet<CloseGenome>();
    }

    /**
     * Add a close genome to one of the sets.
     *
     * @param set			neighbor set
     * @param genomeId		ID of the neighbor genome
     * @param dist			distance to the neighbor genome
     */
    private void addGenome(NavigableSet<CloseGenome> set, String genomeId, double dist) {
        CloseGenome closeness = new CloseGenome(genomeId, dist);
        set.add(closeness);
    }

    /**
     * Record a new genome and its distances.  The genome is added to the new neighbor sets
     * if it is close enough, and the neighbor sets are trimmed if they have gotten too big.
     *
     * @param genomeId		ID of the new genome
     * @param protDist		PheS distance
     * @param rnaDist		16s RNA distance
     */
    public void addGenome(String genomeId, double protDist, double rnaDist) {
        this.addGenome(this.protNeighbors, genomeId, protDist);
        this.addGenome(this.rnaNeighbors, genomeId, rnaDist);
    }

    /**
     * @return a list of the close genomes, sorted by protein distance
     */
    public List<CloseGenome> getProtResults() {
        return new ArrayList<CloseGenome>(this.protNeighbors);
    }

    /**
     * @return a list of the close genomes, sorted by RNA distance
     */
    public List<CloseGenome> getRnaResults() {
        return new ArrayList<CloseGenome>(this.rnaNeighbors);
    }

}
