package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.proteins.kmers.reps.RepGenome;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.MagicMap;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for simple App.
 */
public class AppTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * Test the Role object
     */
    public void testRole() {
        // We need to verify that the normalization works.  We create
        // equivalent roles and insure their checksums are equal.
        String rDesc = "(R)-2-hydroxyacid dehydrogenase, similar to L-sulfolactate dehydrogenase (EC 1.1.1.272)";
        Role rObj1 = new Role("2HydrDehySimiLSulf", rDesc);
        assertEquals("Role ID not stored properly.", "2HydrDehySimiLSulf", rObj1.getId());
        Role rObj2 = new Role("2HydrDehySimiLSulf2", "(R)-2-hydroxyacid dehydrogenase, similar to L-sulfolactate dehydrogenase");
        assertEquals("EC number affects checksum.", rObj1.getChecksum(), rObj2.getChecksum());
        assertEquals("Equal checksums is not equal roles.", rObj1, rObj2);
        assertEquals("Equal checksums does not compare-0 roles.", 0, rObj1.compareTo(rObj2));
        assertEquals("Equal roles have different hashcodes.", rObj1.hashCode(), rObj2.hashCode());
        Role rObj3 = new Role("2HydrDehySimiLSulf3", "(r)-2-hydroxyacid dehydrogenase, similar to L-sulfolactate dehydrogenase");
        assertEquals("Role checksums are case sensitive.", rObj1, rObj3);
        Role rObj4 = new Role("2HydrDehySimiLSulf4", "(R)-2-hydroxyacid dehydrogenase (EC 1.1.1.272), similar to L-sulfolactate dehydrogenase");
        assertEquals("Role checksums affected by EC position.", rObj1, rObj4);
        Role rObj5 = new Role("2HydrDehySimiLSulf5", "(R)-2-hydroxyacid dehydrogenase similar to L-sulfolactate dehydrogenase");
        assertEquals("Role checksum affected by comma.", rObj1, rObj5);
    }

    /**
     * Test the bug with the odd X-heavy protein not having a repeatable similarity.
     * @throws IOException
     */
    public void testXProtein() throws IOException {
        RepGenomeDb testdb = new RepGenomeDb(200);
        File inFile = new File("src/test", "xsmall.fa");
        FastaInputStream inStream = new FastaInputStream(inFile);
        Sequence testSeq = inStream.next();
        RepGenome testGenome = new RepGenome(testSeq);
        testdb.checkGenome(testGenome);
        RepGenomeDb.Representation result = testdb.findClosest(testSeq);
        assertTrue("Should not be an outlier.", result.getSimilarity() > 200);
        double distance = testGenome.distance(testGenome);
        assertEquals("Distance to self is nonzero.", 0.0, distance, 0.0);
        inStream.close();
    }


    /**
     * Test roles IDs.
     * @throws IOException
     */
    public void testRoleMagic() throws IOException {
        File inFile = new File("src/test", "words.txt");
        Scanner roleScanner = new Scanner(inFile);
        roleScanner.useDelimiter("\t|\r\n|\n");
        while (roleScanner.hasNext()) {
            String condensed = roleScanner.next();
            String full = roleScanner.next();
            assertEquals("String did not condense.", condensed, MagicMap.condense(full));
        }
        roleScanner.close();
        // Test registration
        RoleMap magicTable = new RoleMap();
        inFile = new File("src/test", "roles.txt");
        roleScanner = new Scanner(inFile);
        roleScanner.useDelimiter("\t|\r\n|\n");
        while (roleScanner.hasNext()) {
            String roleId = roleScanner.next();
            roleScanner.next();
            assertNull("Wrong ID found", magicTable.get(roleId));
            String roleDesc = roleScanner.next();
            Role newRole = new Role(roleId, roleDesc);
            magicTable.register(newRole);
            assertEquals("Registered ID did not read back.", roleDesc, magicTable.getName(roleId));
        }
        roleScanner.close();
        assertTrue("PheS not found.", magicTable.containsKey("PhenTrnaSyntAlph"));
        assertFalse("Known bad key found.", magicTable.containsKey("PhenTrnaSyntGamm"));
        String modifiedRole = "3-keto-L-gulonate-6-phosphate decarboxylase UlaK putative (L-ascorbate utilization protein D) (EC 4.1.1.85)";
        Role newRole = magicTable.findOrInsert(modifiedRole);
        assertEquals("Wrong ID assigned for modified role.", "3KetoLGulo6PhosDeca6", newRole.getId());
        assertSame("Modified role did not read back.", newRole, magicTable.get("3KetoLGulo6PhosDeca6"));
        modifiedRole = "Unique (new) role string without numbers";
        newRole = magicTable.findOrInsert(modifiedRole);
        assertEquals("Wrong ID assigned for unique role.", "UniqRoleStriWith", newRole.getId());
        assertSame("Unique role did not read back.", newRole, magicTable.get("UniqRoleStriWith"));
        Role findRole = magicTable.findOrInsert(modifiedRole);
        assertSame("Unique role was re-inserted.", newRole, findRole);
        modifiedRole = "Unique (old) role string without numbers";
        findRole = magicTable.findOrInsert(modifiedRole);
        assertTrue("Parenthetical did not change role ID.", findRole != newRole);
        assertEquals("Wrong ID assigned for parenthetical role.", "UniqRoleStriWith2", findRole.getId());
        assertSame("Parenthetical role did not read back.", findRole, magicTable.get("UniqRoleStriWith2"));
        modifiedRole = "Unique (newer) role string without numbers";
        newRole = magicTable.findOrInsert(modifiedRole);
        assertEquals("Wrong ID assigned for newer role.", "UniqRoleStriWith3", newRole.getId());
        assertSame("Parenthetical role did not read back.", newRole, magicTable.get("UniqRoleStriWith3"));
        modifiedRole = "Unique role string 12345 with numbers";
        newRole = magicTable.findOrInsert(modifiedRole);
        assertEquals("Name not stored in role.", modifiedRole, newRole.getName());
        assertEquals("Wrong ID assigned for numbered role.", "UniqRoleStri1234n1", newRole.getId());
        modifiedRole = "Unique role string 12345 with more numbers";
        newRole = magicTable.findOrInsert(modifiedRole);
        assertEquals("Wrong ID assigned for second numbered role.", "UniqRoleStri1234n2", newRole.getId());
        // Test save and load.
        File saveFile = new File("src/test", "roles.ser");
        magicTable.save(saveFile);
        RoleMap newTable = RoleMap.load(saveFile);
        for (Role oldRole : magicTable.values()) {
            newRole = newTable.get(oldRole.getId());
            assertNotNull("Could not find role in loaded table.", newRole);
            assertEquals("Loaded table has wrong role name.", newRole.getName(), oldRole.getName());
            assertEquals("Loaded role has wrong checksum.", newRole, oldRole);
        }
    }

    /**
     * test protein kmers
     */
    public void testKmers() {
        ProteinKmers.setKmerSize(10);
        String myProt1 = "MGMLVPLISKISDLSEEAKACVAACSSVEELDEVRGRYIGRAGALTALLA"; // 50 AA
        String myProt2 = "MDINLFKEELEELAKKAKHMLNETASKNDLEQVKVSLLGKKGLLTLQSAA";
        String myProt3 = "MDINLFKEELKHMLNETASKKGLLTLQSA"; // 30 AA
        ProteinKmers kmer1 = new ProteinKmers(myProt1);
        ProteinKmers kmer2 = new ProteinKmers(myProt2);
        ProteinKmers kmer3 = new ProteinKmers(myProt3);
        assertEquals("Kmer1 has wrong protein.", myProt1, kmer1.getProtein());
        assertEquals("Kmer1 has wrong count.", 41, kmer1.size());
        assertEquals("Kmer1/kmer3 wrong similarity.", 3, kmer2.similarity(kmer3));
        assertEquals("Similarity not commutative.", 3, kmer3.similarity(kmer2));
        assertEquals("Kmer1 too close to kmer2.", 1.0, kmer1.distance(kmer2), 0.0);
        assertEquals("Kmer1 too close to kmer3.", 0.95, kmer2.distance(kmer3), 0.005);
    }

    /**
     * FASTA file test
     * @throws IOException
     */
    public void testFasta() throws IOException {
        File inFasta = new File("src/test", "test.fa");
        FastaInputStream inStream = new FastaInputStream(inFasta);
        ArrayList<Sequence> testSeqs = new ArrayList<Sequence>(5);
        for (Sequence input : inStream) {
            testSeqs.add(input);
        }
        inStream.close();
        assertEquals("Wrong number of sequences.", 5, testSeqs.size());
        Sequence seq = testSeqs.get(0);
        assertEquals("Wrong label for seq 1.", "label1", seq.getLabel());
        assertEquals("Wrong comment for seq 1.", "", seq.getComment());
        assertEquals("Wrong sequence for seq 1.", "tgtgcagcgagccctacagccttggagggaacaacacggactacctgccgctcgtctacccaaagggggtccccctccccaacacaacggttaccagcgtgccgagcg", seq.getSequence());
        seq = testSeqs.get(1);
        assertEquals("Wrong label for seq 2.", "label2", seq.getLabel());
        assertEquals("Wrong comment for seq 2.", "comment2 with spaces", seq.getComment());
        assertEquals("Wrong sequence for seq 2.", "ctcaatgggtccgtagtcggcatcggcagatgtgtataagcagcatgcccgccctctgcag", seq.getSequence());
        seq = testSeqs.get(2);
        assertEquals("Wrong label for seq 3.", "label3", seq.getLabel());
        assertEquals("Wrong comment for seq 3.", "comment3", seq.getComment());
        assertEquals("Wrong sequence for seq 3.", "gtataaagtattggcctgttgag", seq.getSequence());
        seq = testSeqs.get(3);
        assertEquals("Wrong label for seq 4.", "label4", seq.getLabel());
        assertEquals("Wrong comment for seq 4.", "comment4", seq.getComment());
        assertEquals("Wrong sequence for seq 4.", "", seq.getSequence());
        seq = testSeqs.get(4);
        assertEquals("Wrong label for seq 5.", "label5", seq.getLabel());
        assertEquals("Wrong comment for seq 5.", "comment5", seq.getComment());
        assertEquals("Wrong sequence for seq 5.", "ggggccctgaggtcctgagcaagtgggtcggcgagagcgagaaggcgataaggt", seq.getSequence());
        // Write the FASTA back out.
        File outFasta = new File("src/test", "fasta.ser");
        FastaOutputStream outStream = new FastaOutputStream(outFasta);
        outStream.write(testSeqs);
        outStream.close();
        // Read it in and verify.
        inStream = new FastaInputStream(outFasta);
        ArrayList<Sequence> readSeqs = new ArrayList<Sequence>(5);
        while (inStream.hasNext()) {
            readSeqs.add(inStream.next());
        }
        inStream.close();
        assertEquals("Wrong number of records on read back.", 5, readSeqs.size());
        for (int i = 0; i < 5; i++) {
            assertEquals("Compare error at position " + i + " on readback.", testSeqs.get(i), readSeqs.get(i));
        }
    }

    /**
     * Test RepGenome object
     */
    public void testRepGenome() {
        String prot1 = "MSHLAELVASAKAAISQASDVAALDNVRVEYLGKKGHLTLQMTTLRELPPEERPAAGAVI" +
                "NEAKEQVQQALNARKAELESAALNARLAAETIDVSLPGRRIENGGLHPVTRTIDRIESFF" +
                "GELGFTVATGPEIEDDYHNFDALNIPGHHPARADHDTFWFDATRLLRTQTSGVQIRTMKA" +
                "QQPPIRIIAPGRVYRNDYDQTHTPMFHQMEGLIVDTNISFTNLKGTLHDFLRNFFEEDLQ" +
                "IRFRPSYFPFTEPSAEVDVMGKNGKWLEVLGCGMVHPNVLRNVGIDPEVYSGFAFGMGME" +
                "RLTMLRYGVTDLRSFFENDLRFLKQFK";
        RepGenome rep1 = new RepGenome("fig|1005530.3.peg.2208", "Escherichia coli EC4402", prot1);
        assertEquals("Incorrect genome ID parsed.", "1005530.3", rep1.getGenomeId());
        assertEquals("FID not stored.", "fig|1005530.3.peg.2208", rep1.getFid());
        assertEquals("Incorrect name stored.", "Escherichia coli EC4402", rep1.getName());
        assertEquals("Incorrect protein stored.", prot1, rep1.getProtein());
        RepGenome rep2 = new RepGenome("fig|1005530.4.peg.2208", "Escherichia coli EC4402 B", prot1);
        assertFalse("Different genome IDs are still equal.", rep1.equals(rep2));
        assertEquals("Incorrect genome ID in second parse.", "1005530.4", rep2.getGenomeId());
        RepGenome rep3;
        try {
            rep3 = new RepGenome("fig|12345.peg.4", "Invalid genome ID", "");
            fail("Invalid genome ID parsed correctly.");
        } catch (IllegalArgumentException e) {
            // Here the correct exception was thrown.
        }
        rep3 = new RepGenome("fig|1129793.4.peg.2957", "Glaciecola polaris LMG 21857", "MSHLAELVASAKAAISQASDVAALDNVRVEYLGKKGHLTLQMTTLRELPPEERPAAGAVI");
        int sim = rep3.similarity(rep1);
        int sim2 = ((ProteinKmers) rep3).similarity(rep1);
        assertEquals("Similarity depends on subclass used.", sim2, sim);
        // Test genome ordering.
        RepGenome rep4 = new RepGenome("fig|1129793.30.peg.2957", "Test genome 1", "");
        RepGenome rep5 = new RepGenome("fig|129793.30.peg.2957", "Test genome 2", "");
        assertTrue("Rep1 not less than rep2.", rep1.compareTo(rep2) < 0);
        assertTrue("Rep1 not less than rep3.", rep1.compareTo(rep3) < 0);
        assertTrue("Rep3 not less than rep4.", rep3.compareTo(rep4) < 0);
        assertTrue("Rep4 not less than rep5.", rep4.compareTo(rep5) < 0);
        assertTrue("Rep2 not greater than rep1.", rep2.compareTo(rep1) > 0);
        assertTrue("Rep3 not greater than rep1.", rep3.compareTo(rep1) > 0);
        assertTrue("Rep4 not greater than rep3.", rep4.compareTo(rep3) > 0);
        assertTrue("Rep5 not greater than rep4.", rep5.compareTo(rep4) > 0);
        rep3 = new RepGenome("fig|1005530.3.peg.2957", "Glaciecola polaris LMG 21857", "MSHLAELVASAKAAISQASDVAALDNVRVEYLGKKGHLTLQMTTLRELPPEERPAAGAVI");
        assertEquals("Equal genome IDs do not compare 0.", 0, rep1.compareTo(rep3));
    }

    /**
     * test RepGenomeDb
     *
     * @throws IOException
     */
    public void testRepGenomeDb() throws IOException {
        RepGenomeDb repDb = new RepGenomeDb(100);
        assertEquals("Wrong kmer size.", ProteinKmers.kmerSize(), repDb.getKmerSize());
        assertEquals("Wrong key protein.", "Phenylalanyl-tRNA synthetase alpha chain", repDb.getProtName());
        assertEquals("Wrong threshold.", 100, repDb.getThreshold());
        ProteinKmers.setKmerSize(9);
        repDb = new RepGenomeDb(200, "ATP synthase delta chain");
        assertEquals("Wrong kmer size in db2.", 9, repDb.getKmerSize());
        assertEquals("Wrong key protein in db2.", "ATP synthase delta chain", repDb.getProtName());
        assertEquals("Wrong threshold in db2.", 200, repDb.getThreshold());
        // Reset the kmer size.
        ProteinKmers.setKmerSize(10);
        // Process a fasta stream to create a rep-genome DB.
        File fastaFile = new File("src/test", "small.fa");
        FastaInputStream fastaStream = new FastaInputStream(fastaFile);
        repDb = new RepGenomeDb(50);
        repDb.addGenomes(fastaStream);
        fastaStream.close();
        // Find the representative of a genome.
        ProteinKmers testSeq = new ProteinKmers(
                "MSHLAELVASAKAAISQASDVAALDNVRVEYLGKKGHLTLQMTTLRELPPEERPAAGAVI" +
                "NEAKEQVQQALNARKAELESAALNARLAAETIDVSLPGRRIENGGLHPVTRTIDRIESFF" +
                "GELGFTVATGPEIEDDYHNFDALNIPGHHPARADHDTFWFDATRLLRTQTSGVQIRTMKA" +
                "QQPPIRIIAPGRVYRNDYDQTHTPMFHQMEGLIVDTNISFTNLKGTLHDFLRNFFEEDLQ" +
                "IRFRPSYFPFTEPSAEVDVMGKNGKWLEVLGCGMVHPNVLRNVGIDPEVYSGFAFGMGME" +
                "RLTMLRYGVTDLRSFFENDLRFLKQFK");
        RepGenomeDb.Representation result = repDb.findClosest(testSeq);
        assertEquals("E coli not found for E coli protein.", "1005530.3", result.getGenomeId());
        assertTrue("E coli not close enough to E coli protein.", result.getSimilarity() >= 200);
        // Now verify that all the sequences are represented.
        fastaStream = new FastaInputStream(fastaFile);
        for (Sequence inSeq : fastaStream) {
            result = repDb.findClosest(inSeq);
            assertTrue("Genome " + inSeq.getLabel() + " not represented.", result.isRepresented());
        }
        fastaStream.close();
        // Now get all the represented genomes and verify that they are far apart.
        RepGenome[] allReps = repDb.all();
        int n = repDb.size();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                // Verify the similarity thresholds here.
                int compareij = allReps[i].similarity(allReps[j]);
                assertTrue("Genomes " + allReps[i] + " and " + allReps[j] + " are too close.  Score = " + compareij,
                        compareij < repDb.getThreshold());
                // Verify the distance behavior here.
                double distij = allReps[i].distance(allReps[j]);
                double distji = allReps[j].distance(allReps[i]);
                assertEquals("Distance not commutative.", distij, distji, 0.0);
                for (int k = j + 1; k < n; k++) {
                    double distjk = allReps[j].distance(allReps[k]);
                    double distik = allReps[i].distance(allReps[k]);
                    assertTrue("Triangle inequality failure.", distij + distjk >= distik);
                    int compareik = allReps[i].similarity(allReps[k]);
                    if (compareik > compareij && distik > distij) {
                        fail("Greater similarity at greater distance.");
                    } else if (compareik < compareij && distik < distij) {
                        fail("Lesser similarity at lesser distance.");
                    }
                }
            }
        }
        // Save this database.
        File saveFile = new File("src/test", "repdb.ser");
        repDb.save(saveFile);
        // Load it back in.
        ProteinKmers.setKmerSize(6);
        RepGenomeDb newDb = RepGenomeDb.load(saveFile);
        assertEquals("Incorrect protein name loaded.", repDb.getProtName(), newDb.getProtName());
        assertEquals("Incorrect kmer size loaded.", repDb.getKmerSize(), newDb.getKmerSize());
        assertEquals("Incorrect kmer set set.", newDb.getKmerSize(), ProteinKmers.kmerSize());
        assertEquals("Incorrect threshold loaded.", repDb.getThreshold(), newDb.getThreshold());
        assertEquals("Wrong number of genomes loaded.", repDb.size(), newDb.size());
        for (RepGenome oldGenome : allReps) {
            assertEquals("Genome " + oldGenome + " not loaded.", oldGenome, newDb.get(oldGenome.getGenomeId()));
        }

    }
}
