/**
 *
 */
package org.theseed.genome.coupling;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.coupling.FeatureClass.Result;


/**
 * This is the base class for class filters.  A class filter can be used to eliminate classes based on global
 * and/or per-genome criteria.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ClassFilter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ClassFilter.class);

    /**
     * Initialize this class filter.
     *
     * @param processor		parent processor performing the couplings
     *
     * @throws IOException
     */
    public ClassFilter(CouplesProcessor processor) throws IOException { }

    /**
     * Apply the filter to the classes found in a genome.
     *
     * @param results	list of results for a single genome
     *
     * @return the number of class instances removed
     */
    protected abstract int apply(List<FeatureClass.Result> results);

    /**
     * Remove the classes in the specified set from the specified results.
     *
     * @param blacklist		set of IDs for the classes to remove
     * @param results		list of results from which to remove the above classes
     *
     * @return the number of class instances removed
     */
    protected int removeClasses(Set<String> blacklist, List<FeatureClass.Result> results) {
        int retVal = 0;
        for (FeatureClass.Result result : results)
            retVal += result.remove(blacklist);
        log.info("{} class instances removed from {} results.", retVal, results.size());
        return retVal;
    }

    /**
     * @return the displayable name of this filter class
     */
    protected abstract String getName();

    @Override
    public String toString() {
        return this.getName();
    }

    /**
     * Enumeration of filter types.
     */
    public static enum Type {
        NONE, BLACKLIST, LIMITED;

        public ClassFilter create(CouplesProcessor processor) throws IOException {
            ClassFilter retVal = null;
            switch (this) {
            case NONE:
                retVal = new ClassFilter.Null(processor);
                break;
            case BLACKLIST:
                retVal = new BlacklistClassFilter(processor);
                break;
            case LIMITED:
                retVal = new LimitedClassFilter(processor);
                break;
            }
            return retVal;
        }

    }

    /**
     * Dummy class that does not filter.
     */
    public static class Null extends ClassFilter {

        public Null(CouplesProcessor processor) throws IOException {
            super(processor);
        }

        @Override
        protected int apply(List<Result> results) {
            return 0;
        }

        @Override
        protected String getName() {
            return "NONE";
        }


    }

}
