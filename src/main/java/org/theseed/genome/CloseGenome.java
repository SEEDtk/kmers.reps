/**
 *
 */
package org.theseed.genome;

import java.util.Comparator;

/**
 * This object describes a nearby genome.  It sorts with the shortest distance first.
 *
 * @author Bruce Parrello
 */
public class CloseGenome implements Comparable<CloseGenome> {

    // FIELDS
    /** ID of the close genome */
    private String genomeId;
    /** distance to the close genome */
    private double distance;

    /**
     * This is an alternate comparator to sort the furthest distance first.
     */
    public static class Reverse implements Comparator<CloseGenome> {

        @Override
        public int compare(CloseGenome o1, CloseGenome o2) {
            int retVal = Double.compare(o2.getDistance(), o1.getDistance());
            if (retVal == 0)
                retVal = o1.genomeId.compareTo(o2.genomeId);
            return retVal;
        }

    }

    /**
     * Create a close-genome descriptor.
     *
     * @param id		ID of the close genome
     * @param dist		distance to the close genome
     */
    protected CloseGenome(String id, double dist) {
        this.genomeId = id;
        this.distance = dist;
    }

    @Override
    public int compareTo(CloseGenome o) {
        int retVal = Double.compare(this.distance, o.distance);
        if (retVal == 0)
            retVal = this.genomeId.compareTo(o.genomeId);
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.genomeId == null) ? 0 : this.genomeId.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof CloseGenome)) {
            return false;
        }
        CloseGenome other = (CloseGenome) obj;
        if (this.genomeId == null) {
            if (other.genomeId != null) {
                return false;
            }
        } else if (!this.genomeId.equals(other.genomeId)) {
            return false;
        }
        return true;
    }

    /**
     * @return the genome ID
     */
    public String getGenomeId() {
        return this.genomeId;
    }

    /**
     * @return the distance
     */
    public double getDistance() {
        return this.distance;
    }

}
