/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.theseed.genome.CloseGenome;
import org.theseed.genome.CloseSets;

/**
 * This report displays basic statistics for the RNA and seed protein neighborhoods.
 *
 * @author Bruce Parrello
 *
 */
public class StatsNeighborReporter extends NeighborReporter {

    public StatsNeighborReporter(PrintWriter writer) {
        super(writer);
    }

    @Override
    public void openReport() {
        this.println("genome_id\tgenome_name\ttype\tmin\tv100\tq1\tq2\tq3\tmax\tmean\tsdev");
    }

    @Override
    public void recordGenome(String id, String name, CloseSets results) {
        this.output(id, name, "prot", results.getProtResults());
        this.output(id, name, "rna", results.getRnaResults());
    }

    /**
     * Write an output line for a distance list.
     *
     * @param id			genome ID
     * @param name			genome name
     * @param type			type of distance
     * @param results		sorted list of distance results
     */
    private void output(String id, String name, String type, List<CloseGenome> results) {
        DescriptiveStatistics stats = new DescriptiveStatistics(results.stream().mapToDouble(x -> x.getDistance())
                .toArray());
        double q1 = stats.getPercentile(25.0);
        double q2 = stats.getPercentile(50.0);
        double q3 = stats.getPercentile(75.0);
        int i100 = (results.size() <= 100 ? results.size() - 1 : 100);
        double v100 = results.get(i100).getDistance();
        this.print("%s\t%s\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f",
                id, name, type, stats.getMin(), v100, q1, q2, q3, stats.getMax(), stats.getMean(), stats.getStandardDeviation());
    }

    @Override
    public void finish() {
    }

}
