/**
 *
 */
package org.theseed.proteins.kmers.reps;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * @author Bruce Parrello
 *
 */
public class ProtTest {

    /**
     * Test the protein data factory.
     *
     * @throws IOException
     */
    @Test
    public void testProteinFactory() throws IOException {
        // Set up the protein data factory.
        ProteinDataFactory factory = new ProteinDataFactory();
        // Load in the genomes. One of them should fail on genus/species.
        try (TabbedLineReader tabReader = new TabbedLineReader(new File("src/test", "quality.tbl"))) {
            int idCol = tabReader.findField("genome_id");
            int nameCol = tabReader.findField("genome_name");
            int lineageCol = tabReader.findField("taxon_lineage_ids");
            int scoreCol = tabReader.findField("Score");
            for (TabbedLineReader.Line line : tabReader) {
                factory.addGenome(line.get(idCol), line.get(nameCol), line.get(lineageCol), line.getDouble(scoreCol));
            }
            // First pass through.  We expect 1869227.403 to have dropped out because of
            // an insufficient taxonomy.
            Iterator<ProteinData> iter = factory.iterator();
            assertThat(iter.next().getGenomeId(), equalTo("163164.29"));
            assertThat(iter.next().getGenomeId(), equalTo("754252.34"));
            assertThat(iter.next().getGenomeId(), equalTo("813.178"));
            assertThat(iter.next().getGenomeId(), equalTo("813.179"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4260"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4261"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4262"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4263"));
            assertThat(iter.next().getGenomeId(), equalTo("2702.133"));
            assertThat(iter.next().getGenomeId(), equalTo("406327.12"));
            assertThat(iter.next().getGenomeId(), equalTo("571.605"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.5574"));
            assertThat(iter.next().getGenomeId(), equalTo("666.4593"));
            assertThat(iter.next().getGenomeId(), equalTo("1069623.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1313.6345"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.1340"));
            assertThat(iter.next().getGenomeId(), equalTo("1408469.3"));
            assertThat(iter.next().getGenomeId(), equalTo("83555.68"));
            assertThat(iter.next().getGenomeId(), equalTo("2587806.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1733.9145"));
            assertThat(iter.next().getGenomeId(), equalTo("1280.20957"));
            assertFalse(iter.hasNext());
            assertThat(factory.size(), equalTo(21));
            // Now we finalize.  This will drop out 83555.68 due to the lack of
            // a seed protein. We use a small batch size to test the batching.
            factory.finishList(8);
            iter = factory.iterator();
            assertThat(iter.next().getGenomeId(), equalTo("163164.29"));
            assertThat(iter.next().getGenomeId(), equalTo("754252.34"));
            assertThat(iter.next().getGenomeId(), equalTo("813.178"));
            assertThat(iter.next().getGenomeId(), equalTo("813.179"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4260"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4261"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4262"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4263"));
            assertThat(iter.next().getGenomeId(), equalTo("2702.133"));
            assertThat(iter.next().getGenomeId(), equalTo("406327.12"));
            assertThat(iter.next().getGenomeId(), equalTo("571.605"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.5574"));
            assertThat(iter.next().getGenomeId(), equalTo("666.4593"));
            assertThat(iter.next().getGenomeId(), equalTo("1069623.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1313.6345"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.1340"));
            assertThat(iter.next().getGenomeId(), equalTo("1408469.3"));
            assertThat(iter.next().getGenomeId(), equalTo("2587806.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1733.9145"));
            assertThat(iter.next().getGenomeId(), equalTo("1280.20957"));
            assertFalse(iter.hasNext());
            assertThat(factory.size(), equalTo(20));
            // Verify that we have complete information.
            for (ProteinData genomeData : factory) {
                assertNotNull(genomeData.getGenomeName());
                assertNotNull(genomeData.getGenus());
                assertNotNull(genomeData.getSpecies());
                assertNotNull(genomeData.getDna());
                assertNotNull(genomeData.getProtein());
                assertNotNull(genomeData.getFid());
            }
            // Get a specific genome.
            ProteinData genomeData = factory.getGenome("571.605");
            assertThat(genomeData.getSpecies(), equalTo("571"));
            assertThat(genomeData.getGenus(), equalTo("570"));
            assertThat(genomeData.getDna(), equalTo(
                    "atgtcacatctcgcagagctggttgccagtgcgaaggcagccattaacgaggcatcagat" +
                    "gttgctgcgctggacaacgtccgcgtggaatacctgggtaaaaaaggtctcctgaccctt" +
                    "cagatgacgaccctgcgtgagctgcctgctgaagagcgtccggcagccggtgcggttatc" +
                    "aacgaagcgaaagagcaggtccagcaggcgcttaacgcgcgcaagtcagcgctcgaaagc" +
                    "gcagcgctcaacgcgcgtctggcctcggaaaccattgatgtctctctgccggggcgtcgt" +
                    "atcgagaacggtggcctgcatccggtgacccgtaccatcgaccgtattgaaagtttcttc" +
                    "ggtgagctcggttttaccgtcgcgactggcccggagatcgaagatgattatcacaacttc" +
                    "gatgcgctgaatattccaggccaccacccggcacgcgctgaccacgacactttctggttt" +
                    "gatgccacgcgcctgctgcgcacgcaaacatcaggcgtacagatccgcaccatggctaat" +
                    "cagcagccgccaatccgcattattgcccccggccgcgtgtatcgtaacgactacgatcag" +
                    "acgcataccccgatgttccatcagatggaaggtctgatcgttgacactaacatcagcttc" +
                    "accaacctgaagggaacgctgcacgatttcctgcgtaacttctttgaagaagacctgcag" +
                    "attcgttttcgtccgtcctatttcccgttcactgagccgtctgcagaagttgacgtgatg" +
                    "ggtaaaaacggtaaatggctggaagtgctcggctgcggtatggtgcatccaaacgtgctg" +
                    "cgtaacgtgggcatcgatccggaaatctattccggctttgccttcggcatgggtatggag" +
                    "cgcctgaccatgctgcgctatggcgtgaccgacttacgcgcgttcttcgaaaacgatctg" +
                    "cgtttcctcaaacagtttaaataa"));
            assertThat(genomeData.getProtein(), equalTo(
                    "MSHLAELVASAKAAINEASDVAALDNVRVEYLGKKGLLTLQMTTLRELPAEERPAAGAVI" +
                    "NEAKEQVQQALNARKSALESAALNARLASETIDVSLPGRRIENGGLHPVTRTIDRIESFF" +
                    "GELGFTVATGPEIEDDYHNFDALNIPGHHPARADHDTFWFDATRLLRTQTSGVQIRTMAN" +
                    "QQPPIRIIAPGRVYRNDYDQTHTPMFHQMEGLIVDTNISFTNLKGTLHDFLRNFFEEDLQ" +
                    "IRFRPSYFPFTEPSAEVDVMGKNGKWLEVLGCGMVHPNVLRNVGIDPEIYSGFAFGMGME" +
                    "RLTMLRYGVTDLRAFFENDLRFLKQFK"));
            assertThat(genomeData.getGenomeName(), equalTo("Klebsiella oxytoca strain 4928STDY7071345"));
            assertThat(genomeData.getDomain(), equalTo("Bacteria"));
            assertThat(genomeData.getScore(), closeTo(287.622578, 0.00001));
            assertThat(genomeData.getGeneticCode(), equalTo(11));
            // Verify that we find Archaea domains properly.
            genomeData = factory.getGenome("406327.12");
            assertThat(genomeData.getDomain(), equalTo("Archaea"));
            assertThat(genomeData.getGenomeName(), equalTo("Methanococcus vannielii SB"));
            // Verify that we find gc = 4 properly.
            genomeData = factory.getGenome("1408469.3");
            assertThat(genomeData.getGeneticCode(), equalTo(4));
            assertThat(genomeData.getGenus(), equalTo("2093"));
            assertThat(genomeData.getSpecies(), equalTo("36744"));
        }

    }

    /**
     * test the SequenceInfo object
     *
     * @throws IOException
     */
    @Test
    public void testSeqInfo() throws IOException {
        FastaInputStream testStream = new FastaInputStream(new File("src/test", "VapbProt.fa"));
        // Get the first sequence and build the info object.
        Sequence baseSeq = testStream.next();
        SequenceInfo baseInfo = new SequenceInfo(baseSeq);
        assertThat(baseInfo.getMean(), equalTo(0.0));
        assertThat(baseInfo.getSDev(), equalTo(0.0));
        assertThat(baseInfo.getId(), equalTo("fig|1005048.3.peg.2359"));
        // Test the distance to the second sequence and compute the trivial results.
        Sequence newSeq = testStream.next();
        ProteinKmers baseK = new ProteinKmers(baseSeq.getSequence());
        ProteinKmers newK = new ProteinKmers(newSeq.getSequence());
        double dist1 = baseK.distance(newK);
        SequenceInfo newInfo = new SequenceInfo(newSeq);
        baseInfo.storeComparison(newInfo);
        assertThat(baseInfo.getMean(), equalTo(dist1));
        assertThat(baseInfo.getSDev(), equalTo(0.0));
        assertThat(baseInfo.getMin(), equalTo(dist1));
        assertThat(baseInfo.getMax(), equalTo(dist1));
        List<Double> distList = new ArrayList<Double>(350);
        distList.add(dist1);
        for (Sequence testSeq : testStream) {
            newInfo = new SequenceInfo(testSeq);
            dist1 = baseInfo.storeComparison(newInfo);
            distList.add(dist1);
        }
        double mean = baseInfo.getMean();
        double max = baseInfo.getMax();
        double min = baseInfo.getMin();
        double sdev = baseInfo.getSDev();
        testStream.close();
        // Verify that the min, max, and sdev are reasonable.
        int outsideCount = 0;
        int expectedOutside = distList.size() / 4 + 1;
        double sMin = mean - 2 * sdev;
        double sMax = mean + 2 * sdev;
        for (double dist : distList) {
            assertThat(dist, lessThanOrEqualTo(max));
            assertThat(dist, greaterThanOrEqualTo(min));
            if (dist < sMin || dist > sMax) outsideCount++;
        }
        assertThat(outsideCount, lessThan(expectedOutside));
    }


}
