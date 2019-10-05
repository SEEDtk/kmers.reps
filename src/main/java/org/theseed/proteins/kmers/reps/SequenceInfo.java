/**
 *
 */
package org.theseed.proteins.kmers.reps;

import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * This is a utility class that captures information about a sequence for RepMatrixProcessor.
 *
 * It contains the sequence's kmer object plus a summary of the distances to the other
 * sequences.
 *
 * This object can be used in a sort. The comparison is by mean, then standard deviation,
 * then ID.
 *
 * @author Bruce Parrello
 *
 */
public class SequenceInfo implements Comparable<SequenceInfo> {

    // FIELDS

    /** original sequence */
    Sequence original;
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
        this.original = seq;
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
        double retVal = getDistance(other);
        this.addDistance(retVal);
        other.addDistance(retVal);
        return retVal;
    }

    /**
     * @return the distance to another sequence
     *
     * @param other		the sequence whose distance is desired
     */
    public double getDistance(SequenceInfo other) {
        return this.kmers.distance(other.kmers);
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
        double N = this.distCount;
        return (N > 0.0 ? this.distSum / N : 0);
    }

    /**
     * @return the standard deviation of the distance to this sequence
     */
    public double getSDev() {
        double Nm1 = this.distCount - 1;
        double N = this.distCount;
        double retVal = 0.0;
        if (N > 1.0)
            retVal = Math.sqrt(this.distSquareSum / Nm1 - this.distSum * this.distSum / (N * Nm1));
        return retVal;
    }

    /**
     * @return the maximum distance to this sequence
     */
    public double getMax() {
        return this.distMax;
    }

    /**
     * @return the minimum distance to this sequence
     */
    public double getMin() {
        return this.distMin;
    }

    /**
     * @return the ID of this object's sequence
     */
    public String getId() {
        return this.original.getLabel();
    }

    /**
     * @return the original sequence
     */
    public Sequence getSeq() {
        return this.original;
    }

    @Override
    public int compareTo(SequenceInfo o) {
        int retVal = Double.compare(this.distSum, o.distSum);
        if (retVal == 0) {
            retVal = Double.compare(this.getSDev(), o.getSDev());
            if (retVal == 0) {
                retVal = this.getId().compareTo(o.getId());
            }
        }
        return retVal;
    }


}
