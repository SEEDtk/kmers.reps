/**
 *
 */
package org.theseed.genome.coupling;

import java.util.ArrayList;
import java.util.List;

/**
 * This neighbor finder returns all subsequent features within the gap distance.  Note that being on a different
 * strand or contig puts you out of range, and that the incoming result list is sorted by location within strand,
 * so the first result we find beyond the gap distance is a stopping point.
 *
 * @author Bruce Parrello
 *
 */
public class CloseNeighborFinder extends NeighborFinder {

    // FIELDS
    /** initial size for array list */
    private static final int INIT_SIZE = 10;

    public CloseNeighborFinder(BaseCouplingProcessor processor) {
        super(processor);
    }

    @Override
    public List<FeatureClass.Result> getNeighbors(List<FeatureClass.Result> results, int pos) {
        FeatureClass.Result resI = results.get(pos);
        List<FeatureClass.Result> retVal = new ArrayList<FeatureClass.Result>(INIT_SIZE);
        for (int j = pos + 1; j < results.size() && resI.getDistance(results.get(j)) <= this.getMaxGap(); j++)
            retVal.add(results.get(j));
        return retVal;
    }

}
