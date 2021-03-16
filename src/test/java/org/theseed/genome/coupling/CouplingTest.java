/**
 *
 */
package org.theseed.genome.coupling;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.hamcrest.core.IsNull;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;


/**
 * @author Bruce Parrello
 *
 */
public class CouplingTest extends TestCase {

    /**
     * Test feature class results.
     *
     * @throws IOException
     */
    public void testFeatureClasses() throws IOException {
        Genome gto = new Genome(new File("data", "202462.4.gto"));
        Feature feat = gto.getFeature("fig|202462.4.peg.207");
        FeatureClass roleFC = FeatureClass.Type.ROLES.create();
        FeatureClass famFC = FeatureClass.Type.PGFAMS.create();
        FeatureClass randFC = FeatureClass.Type.RANDOM.create();
        FeatureClass.Result res1 = roleFC.getClasses(feat);
        assertThat(res1.getFid(), equalTo(feat.getId()));
        int count = 0;
        for (String role : res1) {
            assertThat(roleFC.getName(role), anyOf(equalTo("SSU ribosomal protein S4p (S9e)"), equalTo("SSU ribosomal protein S4p (S9e), zinc-independent")));
            count++;
        }
        assertThat(count, equalTo(2));
        FeatureClass.Result res3 = famFC.getClasses(feat);
        assertThat(res1, equalTo(res3));
        count = 0;
        for (String fam : res3) {
            assertThat(fam, equalTo("PGF_00049893"));
            count++;
        }
        assertThat(count, equalTo(1));
        FeatureClass.Result res4 = randFC.getClasses(feat);
        count = 0;
        String saved = "";
        for (String fam : res4) {
            boolean found = false;
            for (Feature peg : gto.getPegs()) {
                String pegFam = peg.getPgfam();
                if (pegFam != null && pegFam.contentEquals(fam))
                    found = true;
            }
            assertTrue(found);
            count++;
            saved = fam;
        }
        assertThat(count, equalTo(1));
        res4 = randFC.getClasses(feat);
        count = 0;
        for (String fam : res4) {
            boolean found = false;
            for (Feature peg : gto.getPegs()) {
                String pegFam = peg.getPgfam();
                if (pegFam != null && pegFam.contentEquals(fam))
                    found = true;
            }
            assertTrue(found);
            assertThat(fam, not(equalTo(saved)));
            count++;
        }
        assertThat(count, equalTo(1));
        // Test distances
        feat = gto.getFeature("fig|202462.4.peg.208");
        FeatureClass.Result res2 = famFC.getClasses(feat);
        assertThat(res2, IsNull.notNullValue());
        assertThat(res1.getDistance(res2), equalTo(357));
        assertThat(res2.getDistance(res1), equalTo(357));
        // Overlap distance is less than 0.
        assertThat(res1.getDistance(res3), lessThan(0));
        // test sorting
        assertThat(res1, lessThan(res2));
        // Test different contigs
        feat = new Feature("fig|202462.4.peg.0", feat.getFunction(), "202462.4.con.0003", "+", 190484, 191484);
        res2 = roleFC.getClasses(feat);
        assertThat(res1.getDistance(res2), equalTo(Integer.MAX_VALUE));
        // Test different strands
        feat = gto.getFeature("fig|202462.4.peg.38");
        res2 = famFC.getClasses(feat);
        assertThat(res2, not(equalTo(res1)));
        assertThat(res2.getDistance(res1), equalTo(Integer.MAX_VALUE));
        feat = gto.getFeature("fig|202462.4.peg.569");
        res3 = roleFC.getClasses(feat);
        assertFalse(res3.isGood());
        feat = gto.getFeature("fig|202462.4.rna.12");
        res2 = famFC.getClasses(feat);
        assertFalse(res2.isGood());
    }

    /**
     * Test the ability to read pairs from an output report.
     *
     * @throws IOException
     */
    public void testFeatureClassReads() throws IOException {
        // For families, we rely on the fact that each output file begins with the string representation of a pair.
        FeatureClass fClass = FeatureClass.Type.PGFAMS.create();
        File inFile = new File("data", "coupling10.tbl");
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            for (TabbedLineReader.Line line : inStream) {
                String whole = line.getAll();
                FeatureClass.Pair pairIn = fClass.readPair(line);
                String pairName = pairIn.toString() + "\t";
                assertThat(whole, startsWith(pairName));
            }
        }
        // For roles, the variability of the role IDs may cause the pair to swap.  We have to verify that the pair
        // has matching names.
        fClass = FeatureClass.Type.ROLES.create();
        inFile = new File("data", "coupling10.roles.tbl");
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            for (TabbedLineReader.Line line : inStream) {
                String order1 = line.get(0) + "\t" + line.get(1);
                String order2 = line.get(1) + "\t" + line.get(0);
                FeatureClass.Pair pairIn = fClass.readPair(line);
                assertThat(pairIn.toString(), anyOf(equalTo(order1), equalTo(order2)));
            }
        }
    }

    /**
     * test pairs
     */
    public void testPairs() {
        FeatureClass famFC = FeatureClass.Type.PGFAMS.create();
        FeatureClass.Pair pair1 = famFC.new Pair("x1", "y1");
        FeatureClass.Pair pair2 = famFC.new Pair("y1", "x1");
        FeatureClass.Pair pair3 = famFC.new Pair("x1", "y2");
        assertThat(pair1, equalTo(pair2));
        assertThat(pair2, not(equalTo(pair3)));
        assertThat(pair1.toString(), equalTo("x1\t\ty1\t"));
        assertThat(pair2.toString(), equalTo("x1\t\ty1\t"));
    }

    /**
     * test the result list
     * @throws IOException
     */
    public void testGenomeResult() throws IOException {
        FeatureClass famFC = FeatureClass.Type.PGFAMS.create();
        Genome gto = new Genome(new File("data", "202462.4.gto"));
        // Get the results and insure we got enough.
        List<FeatureClass.Result> results = famFC.getResults(gto);
        assertThat(results.size(), equalTo(gto.getFeatures().size()));
        // Verify the sort.
        for (int i = 0; i < results.size(); i++) {
            FeatureClass.Result resI = results.get(i);
            // Verify the ordering.
            int lastD = Integer.MIN_VALUE;
            for (int j = i + 1; j < results.size(); j++) {
                FeatureClass.Result resJ = results.get(j);
                assertThat(resI.getFid() + " compared to " + resJ.getFid(), resI, lessThan(resJ));
                int thisD = resI.getDistance(resJ);
                assertThat(resI.getFid() + " compared to " + resJ.getFid(), thisD, greaterThanOrEqualTo(lastD));
                lastD = thisD;
            }
            // Verify the feature is valid.
            Feature feat = gto.getFeature(resI.getFid());
            assertThat(feat.getId(), feat, IsNull.notNullValue());
            for (String rClass : resI) {
                assertThat(feat.getId(), rClass, equalTo(feat.getPgfam()));
            }
        }
        // Now we test neighbor finders.
        CouplesProcessor processor = new CouplesProcessor();
        processor.setDefaults();
        NeighborFinder close = NeighborFinder.Type.CLOSE.create(processor);
        NeighborFinder adj = NeighborFinder.Type.ADJACENT.create(processor);
        List<FeatureClass.Result> neighbors = close.getNeighbors(results, 0);
        assertThat(neighbors.size(), equalTo(6));
        assertThat(neighbors.get(0).getFid(), equalTo("fig|202462.4.peg.2"));
        assertThat(neighbors.get(1).getFid(), equalTo("fig|202462.4.peg.4"));
        assertThat(neighbors.get(2).getFid(), equalTo("fig|202462.4.peg.5"));
        assertThat(neighbors.get(3).getFid(), equalTo("fig|202462.4.peg.6"));
        assertThat(neighbors.get(4).getFid(), equalTo("fig|202462.4.peg.7"));
        assertThat(neighbors.get(5).getFid(), equalTo("fig|202462.4.peg.8"));
        neighbors = adj.getNeighbors(results, 0);
        assertThat(neighbors.size(), equalTo(1));
        assertThat(neighbors.get(0).getFid(), equalTo("fig|202462.4.peg.2"));
        neighbors = close.getNeighbors(results, 251);
        assertThat(neighbors.size(), equalTo(0));
        neighbors = adj.getNeighbors(results, 251);
        assertThat(neighbors.size(), equalTo(0));
        neighbors = close.getNeighbors(results, 252);
        assertThat(neighbors.size(), equalTo(2));
        assertThat(neighbors.get(0).getFid(), equalTo("fig|202462.4.peg.607"));
        neighbors = adj.getNeighbors(results, 252);
        assertThat(neighbors.size(), equalTo(1));
        assertThat(neighbors.get(0).getFid(), equalTo("fig|202462.4.peg.607"));
        neighbors = close.getNeighbors(results, 258);
        assertThat(neighbors.size(), equalTo(0));
        neighbors = adj.getNeighbors(results, 258);
        assertThat(neighbors.size(), equalTo(0));


    }

}
