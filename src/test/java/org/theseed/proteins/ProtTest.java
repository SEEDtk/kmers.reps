package org.theseed.proteins;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Iterator;
import org.apache.commons.lang3.StringUtils;
import org.theseed.io.TabbedLineReader;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class ProtTest extends TestCase {

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public ProtTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( ProtTest.class );
    }

    /**
     * Test the role matrix
     */
    public void testRoleMatrix() {
        RoleMatrix testMatrix = new RoleMatrix(2, 5);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(0));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(0));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiHigb"), equalTo(0));
        assertFalse(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        Collection<String> ab = Arrays.asList("AntiHiga", "ToxiHigb");
        Collection<String> ac = Arrays.asList("AntiHiga", "VapcToxiProt");
        Collection<String> abc = Arrays.asList("AntiHiga", "ToxiHigb", "VapcToxiProt");
        Collection<String> abcB = Arrays.asList("AntiHiga", "ToxiHigb", "VapbProt", "VapcToxiProt");
        testMatrix.register("83333.1", ab);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(1));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(0));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        testMatrix.register("83333.1", abc);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(1));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(1));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertTrue(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        testMatrix.register("100226.1", ac);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(2));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(2));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertFalse(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertTrue(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertTrue(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        testMatrix.register("11446.1", ac);
        assertThat(testMatrix.roleCount("AntiHiga"), equalTo(3));
        assertThat(testMatrix.roleCount("ToxiHigb"), equalTo(1));
        assertThat(testMatrix.roleCount("VapbProt"), equalTo(0));
        assertThat(testMatrix.roleCount("VapcToxiProt"), equalTo(3));
        assertTrue(testMatrix.rolePresent("83333.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("100226.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("11446.1", "AntiHiga"));
        assertTrue(testMatrix.rolePresent("83333.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("100226.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("11446.1", "ToxiHigb"));
        assertFalse(testMatrix.rolePresent("83333.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("100226.1", "VapbProt"));
        assertFalse(testMatrix.rolePresent("11446.1", "VapbProt"));
        assertTrue(testMatrix.rolePresent("83333.1", "VapcToxiProt"));
        assertTrue(testMatrix.rolePresent("100226.1", "VapcToxiProt"));
        assertTrue(testMatrix.rolePresent("11446.1", "VapcToxiProt"));
        assertThat(testMatrix.completeness(ab, "83333.1"), equalTo(1.0));
        assertThat(testMatrix.completeness(abc, "83333.1"), equalTo(1.0));
        assertThat(testMatrix.completeness(abcB, "83333.1"), closeTo(0.750, 0.001));
        assertThat(testMatrix.completeness(abcB, "100226.1"), closeTo(0.500, 0.001));
        assertThat(testMatrix.completeness(abcB, "11446.1"), closeTo(0.500, 0.001));
        Collection<String> roles = testMatrix.getMarkerRoles(0.75);
        assertThat(testMatrix.completeness(roles, "83333.1"), greaterThanOrEqualTo(0.75));
        assertThat(testMatrix.completeness(roles, "100226.1"), greaterThanOrEqualTo(0.75));
        assertThat(testMatrix.completeness(roles, "11446.1"), greaterThanOrEqualTo(0.75));
    }

    /**
     * Stress test for role matrix
     * @throws IOException
     */
    public void testRoleMatrixStress() throws IOException {
        // rickettsia file (100 genomes)
        File inFile = new File("src/test", "rickettsia.roles.tbl");
        FileReader fileStream = new FileReader(inFile);
        RoleMatrix stressMatrix = new RoleMatrix(100, 100);
        Set<String> genomes = new HashSet<String>();
        try (BufferedReader reader = new BufferedReader(fileStream)) {
            // Skip the header.
            reader.readLine();
            // Loop through the file.
            for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                String[] fields = StringUtils.splitPreserveAllTokens(line, '\t');
                List<String> roles = Arrays.asList(Arrays.copyOfRange(fields, 2, fields.length));
                stressMatrix.register(fields[0], roles);
                genomes.add(fields[0]);
            }
        }
        Set<String> matrixGenomes = stressMatrix.genomes();
        for (String genome : genomes)
            assertTrue(matrixGenomes.contains(genome));
        assertThat(matrixGenomes.size(), equalTo(genomes.size()));
        Collection<String> universals = stressMatrix.getMarkerRoles(0.90);
        for (String genome : genomes) {
            double completeness = stressMatrix.completeness(universals, genome);
            assertThat(completeness, greaterThanOrEqualTo(0.90));
        }
        universals = stressMatrix.getCommonRoles(0.95);
        for (String role : universals) {
            int frequency = (int) (stressMatrix.roleCount(role) / 0.95);
            assertThat(frequency, greaterThanOrEqualTo(genomes.size()));
        }

    }

    /**
     * Test the protein data factory.
     *
     * @throws IOException
     */
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
            assertThat(iter.next().getGenomeId(), equalTo("813.178"));
            assertThat(iter.next().getGenomeId(), equalTo("813.179"));
            assertThat(iter.next().getGenomeId(), equalTo("163164.29"));
            assertThat(iter.next().getGenomeId(), equalTo("754252.34"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4260"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4261"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4262"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4263"));
            assertThat(iter.next().getGenomeId(), equalTo("2702.133"));
            assertThat(iter.next().getGenomeId(), equalTo("406327.12"));
            assertThat(iter.next().getGenomeId(), equalTo("571.605"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.5574"));
            assertThat(iter.next().getGenomeId(), equalTo("666.4593"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.1340"));
            assertThat(iter.next().getGenomeId(), equalTo("1069623.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1313.6345"));
            assertThat(iter.next().getGenomeId(), equalTo("1408469.3"));
            assertThat(iter.next().getGenomeId(), equalTo("2587806.3"));
            assertThat(iter.next().getGenomeId(), equalTo("83555.68"));
            assertThat(iter.next().getGenomeId(), equalTo("1733.9145"));
            assertThat(iter.next().getGenomeId(), equalTo("1280.20957"));
            assertFalse(iter.hasNext());
            assertThat(factory.size(), equalTo(21));
            // Now we finalize.  This will drop out 83555.68 due to the lack of
            // a seed protein. We use a small batch size to test the batching.
            factory.finishList(8);
            iter = factory.iterator();
            assertThat(iter.next().getGenomeId(), equalTo("813.178"));
            assertThat(iter.next().getGenomeId(), equalTo("813.179"));
            assertThat(iter.next().getGenomeId(), equalTo("163164.29"));
            assertThat(iter.next().getGenomeId(), equalTo("754252.34"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4260"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4261"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4262"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.4263"));
            assertThat(iter.next().getGenomeId(), equalTo("2702.133"));
            assertThat(iter.next().getGenomeId(), equalTo("406327.12"));
            assertThat(iter.next().getGenomeId(), equalTo("571.605"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.5574"));
            assertThat(iter.next().getGenomeId(), equalTo("666.4593"));
            assertThat(iter.next().getGenomeId(), equalTo("1639.1340"));
            assertThat(iter.next().getGenomeId(), equalTo("1069623.3"));
            assertThat(iter.next().getGenomeId(), equalTo("1313.6345"));
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


}
