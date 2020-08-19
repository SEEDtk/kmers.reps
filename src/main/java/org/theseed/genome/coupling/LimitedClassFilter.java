/**
 *
 */
package org.theseed.genome.coupling;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.theseed.counters.CountMap;

/**
 * This filter removes classes that occur too often in the genome.  No processing is done during construction.
 * Instead, a filter set is built on a per-genome basis.
 *
 * @author Bruce Parrello
 *
 */
public class LimitedClassFilter extends ClassFilter {

    // FIELDS
    /** class limit */
    private int classLimit;

    public LimitedClassFilter(CouplesProcessor processor) throws IOException {
        super(processor);
        this.classLimit = processor.getClassLimit();
    }

    @Override
    protected int apply(List<FeatureClass.Result> results) {
        // Parse the result set and count the classes.
        CountMap<String> classCounts = new CountMap<String>();
        for (FeatureClass.Result result : results) {
            for (String classId : result)
                classCounts.count(classId);
        }
        // Pull out the classes that occur more than the limit level.
        Set<String> blackList = classCounts.counts().stream().filter(k -> (k.getCount() > this.classLimit))
                .map(k -> k.getKey()).collect(Collectors.toSet());
        // Remove them from the results.
        return this.removeClasses(blackList, results);
    }

    @Override
    protected String getName() {
        return String.format("LIMITED to classes having no more than %d occurrences per genome", this.classLimit);
    }

}
