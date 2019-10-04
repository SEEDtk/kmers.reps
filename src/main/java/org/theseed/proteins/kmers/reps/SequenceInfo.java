/**
 *
 */
package org.theseed.proteins.kmers.reps;

import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.sequence.Sequence;

// TODO write tests for this

/**
 * This is a utility class that captures information about a sequence for RepMatrixProcessor.
 *
 * It contains the sequence's kmer object plus a summary of the distances to the other
 * sequences.
 *
 * @author Bruce Parrello
 *
 */
public class SequenceInfo {

    // FIELDS

    /** ID of this sequence */
    String	seqId;
    /** protein kmer object */
    ProteinKmers kmers;
    /** number of sequences compared */
    int distCount;
    /** sum of the distances */
    double distSum;
    /** sum of the squares of the distances */
    double distSquareSum;
    /** maximum distance encountered */
    double distMax;
    /** minimum distance encountered */
    double distMin;

    /**
     * Create a sequence info object for a particular sequence.
     *
     * @param seq		sequence for which the info should be initialized
     */
    public SequenceInfo(Sequence seq) {
        this.seqId = seq.getLabel();
        this.kmers = new ProteinKmers(seq.getSequence());
        this.distCount = 0;
        this.distSum = 0;
        this.distSquareSum = 0;
        this.distMax = 0;
        this.distMin = Double.MAX_VALUE;
    }

    /**
     * Compare this sequence with another and update all the relevant counts.
     *
     * @param other		sequence to compare to this one
     *
     * @return the distance between the sequences
     */
    public double storeComparison(SequenceInfo other) {
        double retVal = this.kmers.distance(other.kmers);
        this.addDistance(retVal);
        other.addDistance(retVal);
        return retVal;
    }

    /**
     * Record a new distance in this object.
     *
     * @param dist		distance to another sequence
     */
    private void addDistance(double dist) {
        this.distCount++;
        this.distSum += dist;
        this.distSquareSum += dist * dist;
        if (dist > this.distMax) this.distMax = dist;
        if (dist < this.distMin) this.distMin = dist;
    }

    /**
     * @return the mean distance to this sequence
     */
    public double getMean() {
        return this.distSum / this.distCount;
    }

    /**
     * @return the standard deviation of the distance to this sequence
     */
    public double getSDev() {
        double Nm1 = this.distCount - 1;
        double N = this.distCount;
        double retVal = Math.sqrt(this.distSquareSum / Nm1 - this.distSum * this.distSum / (N * Nm1));
        return retVal;
    }

}
