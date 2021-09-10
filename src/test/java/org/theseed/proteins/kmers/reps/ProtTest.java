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
        try (TabbedLineReader tabReader = new TabbedLineReader(new File("data", "quality.tbl"))) {
            int idCol = tabReader.findField("Genome");
            int nameCol = tabReader.findField("Name");
            int lineageCol = tabReader.findField("Taxonomy");
            int scoreCol = tabReader.findField("Score");
            for (TabbedLineReader.Line line : tabReader) {
                factory.addGenome(line.get(idCol), line.get(nameCol), line.get(lineageCol), line.getDouble(scoreCol));
            }
            // First pass through.  We expect 100480.3 to have dropped out because of
            // an insufficient taxonomy.
            Iterator<ProteinData> iter = factory.iterator();
            assertThat(iter.next().getGenomeId(), equalTo("99287.12"));
            assertThat(iter.next().getGenomeId(), equalTo("309807.25"));
            assertThat(iter.next().getGenomeId(), equalTo("932213.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1001994.6"));
            assertThat(iter.next().getGenomeId(), equalTo("1408469.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1003189.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1003185.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1003186.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1003187.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1004951.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1003188.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1003191.3"));
            assertThat(iter.next().getGenomeId(), equalTo("571.605"));
            assertThat(iter.next().getGenomeId(), equalTo("1003192.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1004952.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1955272.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1004901.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1004836.5"));
            assertThat(iter.next().getGenomeId(), equalTo("37372.4"));
            assertFalse(iter.hasNext());
            assertThat(factory.size(), equalTo(19));
            // Now we finalize.  This will also reorder 37372.4 and 1955272.3 due to
            // bad SSU rRNA.  We use a small batch size to test the batching.
            factory.finishList(8);
            iter = factory.iterator();
            ProteinData.Rating oldRating = ProteinData.Rating.NCBI_REF;
            double oldScore = Double.POSITIVE_INFINITY;
            while (iter.hasNext()) {
                ProteinData genome = iter.next();
                String label = genome.getGenomeName();
                ProteinData.Rating rating = genome.getRating();
                double score = genome.getScore();
                int ratingComp = oldRating.compareTo(rating);
                assertThat(label, ratingComp, lessThanOrEqualTo(0));
                if (ratingComp == 0)
                    assertThat(label, oldScore, greaterThanOrEqualTo(score));
                oldRating = rating;
                oldScore = score;
            }
            assertThat(factory.size(), equalTo(19));
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
            assertThat(genomeData.getSsuSequence(), equalTo(
                    "cttaaattgaagagtttgatcatggctcagattgaacgctggcggcaggcctaacacatg"
                    + "caagtcgaacggtagcacagagagcttgctctcgggtgacgagtggcggacgggtgagta"
                    + "atgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgca"
                    + "taacgtcgcaagaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagat"
                    + "gggattagcttgtaggtgaggtaacggctcacctaggcgacgatccctagctggtctgag"
                    + "aggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtg"
                    + "gggaatattgcacaatgggcgcaagcctgatgcagccatgccgcgtgtatgaagaaggcc"
                    + "ttcgggttgtaaagtactttcagcggggaggaagggagtgaggttaataacctcattcat"
                    + "tgacgttacccgcagaagaagcaccggctaactccgtgccagcagccgcggtaatacgga"
                    + "gggtgcaagcgttaatcggaattactgggcgtaaagcgcacgcaggcggtctgtcaagtc"
                    + "ggatgtgaaatccccgggctcaacctgggaactgcattcgaaactggcaggctggagtct"
                    + "tgtagaggggggtagaattccaggtgtagcggtgaaatgcgtagagatctggaggaatac"
                    + "cggtggcgaaggcggccccctggacaaagactgacgctcaggtgcgaaagcgtggggagc"
                    + "aaacaggattagataccctggtagtccacgctgtaaacgatgtcgacttggaggttgttc"
                    + "ccttgaggagtggcttccggagctaacgcgttaagtcgaccgcctggggagtacggccgc"
                    + "aaggttaaaactcaaatgaattgacgggggcccgcacaagcggtggagcatgtggtttaa"
                    + "ttcgatgcaacgcgaagaaccttacctactcttgacatccanccgcctggggagtacggc"
                    + "cgcaaggttaaaactcaaatgaattgacgggggcccgcacaagcggtggagcatgtggtt"
                    + "taattcgatgcaacgcgaagaaccttacctactcttgacatccacggaatttggcagaga"
                    + "tgccttagtgccttcgggaaccgtgagacaggtgctgcatggctgtcgtcagctcgtgtt"
                    + "gtgaaatgttgggttaagtcccgcaacgagcgcaacccttatcctttgttgccagcggtc"
                    + "cggccgggaactcaaaggagactgccagtgataaactggaggaaggtggggatgacgtca"
                    + "agtcatcatggcccttacgagtagggctacacacgtgctacaatggcatatacaaagaga"
                    + "agcgacctcgcgagagcaagcggacctcataaagtatgtcgtagtccggattggagtctg"
                    + "caactcgactccatgaagtcggaatcgctagtaatcgtggatcagaatgccacggtgaat"
                    + "acgttcccgggccttgtacacaccgcccgtcacaccatgggagtgggttgcaaaagaagt"
                    + "aggtagcttaaccttcgggagggcgcttaccactttgtgattcatgactggggtgaagtc"
                    + "gtaacaaggtaaccgtaggggaacctgcggttggatcacctccttacc"));
            assertThat(genomeData.getGenomeName(), equalTo("Klebsiella oxytoca strain 4928STDY7071345"));
            assertThat(genomeData.getDomain(), equalTo("Bacteria"));
            assertThat(genomeData.getScore(), closeTo(282.73, 0.01));
            assertThat(genomeData.getGeneticCode(), equalTo(11));
            // Verify that we find Archaea domains properly.
            genomeData = factory.getGenome("1001994.6");
            assertThat(genomeData.getDomain(), equalTo("Archaea"));
            assertThat(genomeData.getGenomeName(), equalTo("Candidatus Nitrosoarchaeum koreensis MY1"));
            // Verify that we find gc = 4 properly.
            genomeData = factory.getGenome("1408469.3");
            assertThat(genomeData.getGeneticCode(), equalTo(4));
            assertThat(genomeData.getGenus(), equalTo("2093"));
            assertThat(genomeData.getSpecies(), equalTo("36744"));
            // Verify that the bad guys are bad.
            assertThat(factory.getGenome("1955272.3").getRating(), equalTo(ProteinData.Rating.BAD_SSU));
            assertThat(factory.getGenome("37372.4").getRating(), equalTo(ProteinData.Rating.BAD_SSU));
        }

    }

    /**
     * test the SequenceInfo object
     *
     * @throws IOException
     */
    @Test
    public void testSeqInfo() throws IOException {
        FastaInputStream testStream = new FastaInputStream(new File("data", "VapbProt.fa"));
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
