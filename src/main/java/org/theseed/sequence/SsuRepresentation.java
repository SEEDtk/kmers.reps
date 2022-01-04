/**
 *
 */
package org.theseed.sequence;

/**
 * This class is used to return results from an SSU representation scan.  It describes the
 * closest genome, its distance, and its similarity count.
 *
 * @author Bruce Parrello
 *
 */
public class SsuRepresentation {

    // FIELDS
    /** ID of the close genome */
    private String genomeId;
    /** name of the close genome */
    private String genomeName;
    /** similarity of the close genome */
    private int hitCount;
    /** distance to the close genome */
    private double distance;

    /**
     * Create a blank, empty SSU representation.
     */
    public SsuRepresentation() {
        this.genomeId = "";
        this.genomeName = "";
        this.hitCount = 0;
        this.distance = 1.0;
    }

    /**
     * Update this representation if a new genome is better.
     *
     * @param input			input genome
     * @param candidate		potential representative genome
     *
     * @return TRUE if the new genome is better, else FALSE
     */
    public boolean checkSsu(SsuGenome input, SsuGenome candidate) {
        boolean retVal = false;
        // Compute the similarity and distance.
        int sim = input.similarity(candidate);
        double dist = SequenceKmers.distance(sim, input, candidate);
        // If the candidate is closer, save its information.
        if (dist < this.distance) {
            this.genomeId = candidate.getGenomeId();
            this.genomeName = candidate.getName();
            this.hitCount = sim;
            this.distance = dist;
            retVal = true;
        }
        return retVal;
    }

    /**
     * @return the ID of the close genome
     */
    public String getGenomeId() {
        return this.genomeId;
    }

    /**
     * @return the name of the close genome
     */
    public String getGenomeName() {
        return this.genomeName;
    }

    /**
     * @return the similarity score of the close genome
     */
    public int getHitCount() {
        return this.hitCount;
    }

    /**
     * @return the distance of the close genome
     */
    public double getDistance() {
        return this.distance;
    }

    /**
     * @return TRUE if there is no genome found in this representation
     */
    public boolean isEmpty() {
        return this.genomeId.isEmpty();
    }

}
