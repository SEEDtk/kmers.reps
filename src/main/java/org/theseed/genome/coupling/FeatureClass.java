/**
 *
 */
package org.theseed.genome.coupling;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;

/**
 * This object is used to classify features for functional coupling.  It takes as input a genome feature and
 * produces as output a list of class ID strings.  So, if we are coupling on protein families, it will
 * return a list of the feature's protein families.  If we are coupling on roles, it will return a list
 * of the feature's roles.  An empty list indicates the feature is not valid for consideration.  Note
 * that we expect the number of classes to be small.  There will usually be one, sometimes two or three.
 *
 * @author Bruce Parrello
 *
 */
public abstract class FeatureClass {

    /**
     * @return a (possibly empty) list of classifications for a feature
     *
     * @param feat	feature of interest
     */
    public abstract Result getClasses(Feature feat);

    /**
     * @return the name of a specified class to appear on reports
     *
     * @param classId	ID of the class of interest
     */
    public abstract String getName(String classId);

    /**
     * @return the headers for a pair of classes in the output report
     */
    public abstract String getHeadings();

    /**
     * Insure that names are known for the specified class IDs.
     */
    public abstract void cacheNames(Collection<String> classes);

    /**
     * @return the pair of classes represented by the names in the specified report line
     *
     * @param line	output line from a coupling report
     */
    public abstract Pair readPair(TabbedLineReader.Line line);

    /**
     * @return a sorted list of the classifications for the features in a genome
     *
     * @param genome	genome of interest
     */
    public List<Result> getResults(Genome genome) {
        Collection<Feature> feats = genome.getFeatures();
        List<Result> retVal = new ArrayList<Result>(feats.size());
        for (Feature feat : feats) {
            Result res = this.getClasses(feat);
            retVal.add(res);
        }
        // Sort the results by location.
        retVal.sort(null);
        return retVal;
    }

    /**
     * This class represents a list of feature classes.  It serves as the result object when analyzing features.
     */
    public static class Result implements Comparable<Result>, Iterable<String> {

        // FIELDS
        /** ID of the feature that produced the classifications */
        private String fid;
        /** location of the feature */
        private Location loc;
        /** list of classes found */
        private Set<String> classes;

        /**
         * Construct an empty result for a specified feature.
         */
        protected Result(Feature feat) {
            this.fid = feat.getId();
            this.loc = feat.getLocation();
            this.classes = new HashSet<String>(5);
        }

        /**
         * Add a class to this result set.
         *
         * @param classId	ID of the class to add
         */
        protected void add(String classId) {
            this.classes.add(classId);
        }

        /**
         * @return the distance between the features represented by two results.
         *
         * @param other		other result to check
         */
        public int getDistance(Result other) {
            int retVal;
            if (this.loc.getDir() != other.loc.getDir())
                retVal = Integer.MAX_VALUE;
            else
                retVal = this.loc.distance(other.loc);
            return retVal;
        }

        /**
         * We want to sort results based on location, with strand coming before begin point.
         */
        @Override
        public int compareTo(Result o) {
            int retVal = this.loc.getContigId().compareTo(o.loc.getContigId());
            if (retVal == 0) {
                retVal = this.loc.getDir() - o.loc.getDir();
                if (retVal == 0) {
                    retVal = this.loc.getBegin() - o.loc.getBegin();
                    if (retVal == 0) {
                        retVal = this.fid.compareTo(o.fid);
                    }
                }
            }
            return retVal;
        }

        /**
         * The ID of this object is its feature ID.
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((fid == null) ? 0 : fid.hashCode());
            return result;
        }

        /**
         * The ID of this object is its feature ID.
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Result)) {
                return false;
            }
            Result other = (Result) obj;
            if (fid == null) {
                if (other.fid != null) {
                    return false;
                }
            } else if (!fid.equals(other.fid)) {
                return false;
            }
            return true;
        }

        /**
         * @return an iterator for all the classes in this result
         */
        @Override
        public Iterator<String> iterator() {
            return this.classes.iterator();
        }

        /**
         * @return TRUE if the result has any classes in it
         */
        public boolean isGood() {
            return this.classes.size() > 0;
        }

        /**
         * @return the feature ID
         */
        public String getFid() {
            return fid;
        }

        /**
         * Remove all classes in the specified set from our class list.
         *
         * @param blacklist		set of classes to remove
         *
         * @return the number of classes removed
         */
        public int remove(Set<String> blacklist) {
            int retVal = 0;
            Iterator<String> iter = this.classes.iterator();
            while (iter.hasNext()) {
                String curr = iter.next();
                if (blacklist.contains(curr)) {
                    iter.remove();
                    retVal++;
                }
            }
            return retVal;
        }

    }

    /**
     * This represents an unordered pair of feature classifications.
     *
     * The name of a pair is the names of the two classes, tab-delimited.
     */
    public class Pair {

        // FIELDS
        /** first class */
        private String class1;
        /** second class */
        private String class2;

        /**
         * Construct a pair.
         */
        public Pair(String classA, String classB) {
            // We do a fast sort to insure the order doesn't matter.
            if (classA.compareTo(classB) < 0) {
                this.class1 = classA;
                this.class2 = classB;
            } else {
                this.class1 = classB;
                this.class2 = classA;
            }
        }

        /**
         * @return the name of this pair for reports
         */
        public String toString() {
            return FeatureClass.this.getName(this.class1) + "\t" + FeatureClass.this.getName(this.class2);
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((class1 == null) ? 0 : class1.hashCode());
            result = prime * result + ((class2 == null) ? 0 : class2.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Pair)) {
                return false;
            }
            Pair other = (Pair) obj;
            if (class1 == null) {
                if (other.class1 != null) {
                    return false;
                }
            } else if (!class1.equals(other.class1)) {
                return false;
            }
            if (class2 == null) {
                if (other.class2 != null) {
                    return false;
                }
            } else if (!class2.equals(other.class2)) {
                return false;
            }
            return true;
        }

        /**
         * @return the first class in the pair
         */
        public String getClass1() {
            return this.class1;
        }

        /**
         * @return the second class in the pair
         */
        public String getClass2() {
            return this.class2;
        }

    }

    /**
     * This enumerates the valid classification types.
     */
    public static enum Type {
        ROLES, PGFAMS, RANDOM;

        /**
         * Create the classifier for this classification type.
         */
        public FeatureClass create() {
            FeatureClass retVal = null;
            switch (this) {
            case ROLES :
                retVal = new RoleFeatureClass();
                break;
            case PGFAMS :
                retVal = new FamilyFeatureClass();
                break;
            case RANDOM :
                retVal = new RandomFeatureClass();
                break;
            }
            return retVal;
        }
    }
}
