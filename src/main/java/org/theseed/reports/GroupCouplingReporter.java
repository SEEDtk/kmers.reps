/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Genome;
import org.theseed.genome.coupling.BaseCouplingProcessor;
import org.theseed.genome.coupling.FeatureClass;

/**
 * This is the simplest possible sort of coupling report.  For each output class pair, we simply list the genomes found.
 *
 * @author Bruce Parrello
 */
public class GroupCouplingReporter extends CouplingReporter {

    /**
     * Construct a report to a specified output stream.
     *
     * @param output		target output stream
     * @param processor		parent coupling processor
     */
    public GroupCouplingReporter(OutputStream output, BaseCouplingProcessor processor) {
        super(output, processor);
    }

    @Override
    public void register(Genome genome) { }

    @Override
    protected String getScoreHeadings() {
        return "genomes";
    }

    @Override
    public void writePairLine(FeatureClass.Pair pair, FeatureClass.PairData genomeData) {
        this.print("%s\t%s", pair.toString(), StringUtils.join(genomeData.getGenomes(), ','));
    }

}
