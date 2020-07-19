/**
 *
 */
package org.theseed.genome.coupling;

import java.util.Collections;
import java.util.List;

import org.theseed.genome.coupling.FeatureClass;

/**
 * This is the simplest neighbor finder.  It returns the next feature if it is within the gap distance,
 * and an empty set otherwise.
 * @author Bruce Parrello
 *
 */
public class AdjacentNeighborFinder extends NeighborFinder {

    public AdjacentNeighborFinder(BaseCouplingProcessor processor) {
        super(processor);
    }

    @Override
    public List<FeatureClass.Result> getNeighbors(List<FeatureClass.Result> results, int pos) {
        FeatureClass.Result resJ = results.get(pos + 1);
        List<FeatureClass.Result> retVal;
        if (resJ == null || resJ.getDistance(results.get(pos)) > this.getMaxGap())
            retVal = Collections.emptyList();
        else
            retVal = Collections.singletonList(resJ);
        return retVal;
    }

}
