/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * This class loads a list of sequences into a named sequence collection.  The sequences can then be used to
 * determine how far a sequence is from a specific collection.
 *
 * @author Bruce Parrello
 *
 */
public class KmerCollectionGroup {

    // FIELDS
    /** map of list names to sequence lists */
    HashMap<String, Collection<ProteinKmers>> sequenceLists;

    /**
     * Create new, blank sequence-list collection
     */
    public KmerCollectionGroup() {
        this.sequenceLists = new HashMap<>();
    }

    /**
     * Add a sequence to one of the lists.
     *
     * @param seq		sequence to add
     * @param label		name to give to the sequence's list
     */
    public void addSequence(Sequence seq, String label) {
        // Create the protein kmers object.
        ProteinKmers kmers = new ProteinKmers(seq.getSequence());
        // Find the list for this label.
        Collection<ProteinKmers> target = this.sequenceLists.get(label);
        if (target == null ) {
            target = new ArrayList<ProteinKmers>();
            this.sequenceLists.put(label, target);
        }
        // Add this sequence to the list.
        target.add(kmers);
    }

    /**
     * Find the distance from a sequence to a label's collection.
     *
     * @param seq		sequence to test
     * @param label		list against which to test it
     */
    public double getDistance(Sequence seq, String label) {
        // Start with the max possible distance.
        double retVal = 1.0;
        // Find the list for this label.
        Collection<ProteinKmers> target = this.sequenceLists.get(label);
        if (target != null) {
            // Create kmers for the sequence.
            ProteinKmers kmers = new ProteinKmers(seq.getSequence());
            // Search for the best distance.
            for (ProteinKmers protein : target) {
                double dist = protein.distance(kmers);
                if (dist < retVal) retVal = dist;
            }
        }
        // Return the best distance.
        return retVal;
    }

    /**
     * Result class for returning the best group.
     */
    public static class Result {
        private double distance;
        private String group;

        private Result() {
            this.distance = 1.0;
            this.group = null;
        }

        /**
         * @return the best distance
         */
        public double getDistance() {
            return distance;
        }

        /**
         * @return the best group
         */
        public String getGroup() {
            return group;
        }

        /**
         * Update this object if we have a better result.
         *
         * @param grp	relevant group
         * @param dist	new distance
         */
        protected void merge(String grp, double dist) {
            if (dist < this.distance) {
                this.group = grp;
                this.distance = dist;
            }
        }

    }

    /**
     * Find the group closest to the specified sequence and return the name and distance.
     *
     * @param seq	sequence to check
     *
     * @return a Result object containing the group name and distance
     */
    public Result getBest(Sequence seq) {
        Result retVal = new Result();
        for (String grp : this.sequenceLists.keySet()) {
            double distance = this.getDistance(seq, grp);
            retVal.merge(grp, distance);
        }
        return retVal;
    }

}
