/**
 *
 */
package org.theseed.genome;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.theseed.test.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
import org.theseed.utils.ParseFailureException;

/**
 * Test genome descriptors for RepGen
 *
 * @author Bruce Parrello
 *
 */
public class GenomeDescriptorTest {

    @Test
    public void testDescriptor() throws IOException, ParseFailureException {
        Genome g1 = new Genome(new File("data", "1280.25101.gto"));
        Genome g2 = new Genome(new File("data", "1313.16015.gto"));
        Genome g4 = new Genome(new File("data", "1280.21158.gto"));
        GenomeDescriptor g1d = new GenomeDescriptor(g1);
        GenomeDescriptor g2d = new GenomeDescriptor(g2);
        GenomeDescriptor g4d = new GenomeDescriptor(g4);
        GenomeDescriptor g4da = new GenomeDescriptor(g4.getId(), g4.getName(),
                "MSEQQTMSELKQQALVDINEANDERALQEVKVKYLGKKGSVSGLMKLMKDLPNEDKPAFGQKVNELRQTIQNELDERQQMLVKEKLNKQLAEETIDVSLPGRHIEIGSKHPLTRTIEEIEDLFLGLGYEIVNGYEVEQDHYNFEMLNLPKSHPARDMQDSFYITDEILLRTHTSPVQARTMESRHGQGPVKIICPGKVYRRDSDDATHSHQFTQIEGLVVDKNVKMSDLKGTLELLAKKLFGADREIRLRPSYFPFTEPSVEVDVSCFKCKGKGCNVCKHTGWIEILGAGMVHPNVLEMAGFDSSEYSGFAFGMGPDRIAMLKYGIEDIRHFYTNDVRFLDQFKAVEDRGDM",
                g4.getSsuRRna());
        assertThat(g1d.getId(), equalTo("1280.25101"));
        assertThat(g1d.getName(), equalTo("Staphylococcus aureus strain 12593_2_37"));
        assertThat(g4da.getId(), equalTo(g4.getId()));
        assertThat(g4da.getName(), equalTo(g4.getName()));
        assertThat(g4d.getRnaDistance(g4da), equalTo(0.0));
        assertThat(g4da.getRnaDistance(g4d), equalTo(0.0));
        assertThat(g4d.getSeedDistance(g4da), equalTo(0.0));
        assertThat(g4da.getSeedDistance(g4d), equalTo(0.0));
        int g1_g2_seedSim = g1d.getSeedSim(g2d);
        assertThat(g2d.getSeedSim(g1d), equalTo(g1_g2_seedSim));
        int g1_g4_seedSim = g1d.getSeedSim(g4d);
        assertThat(g1_g2_seedSim, lessThan(g1_g4_seedSim));
        int g1_g2_rnaSim = g1d.getRnaSim(g2d);
        assertThat(g2d.getRnaSim(g1d), equalTo(g1_g2_rnaSim));
        int g1_g4_rnaSim = g1d.getRnaSim(g4d);
        assertThat(g1_g2_rnaSim, lessThan(g1_g4_rnaSim));
        double g1_g2_seedDistance = g1d.getSeedDistance(g2d);
        assertThat(g2d.getSeedDistance(g1d), equalTo(g1_g2_seedDistance));
        double g1_g4_seedDistance = g1d.getSeedDistance(g4d);
        assertThat(g1_g2_seedDistance, greaterThan(g1_g4_seedDistance));
        double g1_g2_rnaDistance = g1d.getRnaDistance(g2d);
        assertThat(g2d.getRnaDistance(g1d), equalTo(g1_g2_rnaDistance));
        double g1_g4_rnaDistance = g1d.getRnaDistance(g4d);
        assertThat(g1_g2_rnaDistance, greaterThan(g1_g4_rnaDistance));
    }

    public static final String[] GENOME_IDS = new String[] { "103621.4", "1036677.3", "1036778.3", "1280.21158", "1313.16015" };

    @Test
    public void testSet() throws IOException, ParseFailureException {
        GenomeDescriptorSet refSet = new GenomeDescriptorSet();
        for (String genomeId : GENOME_IDS) {
            File genomeFile = new File("data", genomeId + ".gto");
            Genome genome = new Genome(genomeFile);
            refSet.add(genome);
        }
        Genome testGenome = new Genome(new File("data", "1280.25101.gto"));
        GenomeDescriptor testDescriptor = new GenomeDescriptor(testGenome);
        List<GenomeDescriptorSet.CloseGenome> resultList = new ArrayList<>(6);
        for (GenomeDescriptorSet.FinderType type : GenomeDescriptorSet.FinderType.values()) {
            GenomeDescriptorSet.CloseGenome result = refSet.findClosest(testDescriptor, type);
            assertThat(result.getGenomeId(), equalTo("1280.21158"));
            assertThat(result.getGenomeName(), equalTo("Staphylococcus aureus strain BPH2869"));
            assertThat(result.getProximity(), greaterThan(0.0));
            resultList.add(result);
        }
        assertThat(GenomeDescriptorSet.CloseGenome.test(resultList, new int[] { 0, 1, 2, 3} ), isTrue());
        resultList.add(new GenomeDescriptorSet.CloseGenome());
        assertThat(GenomeDescriptorSet.CloseGenome.test(resultList, new int[] { 4, 0, 1, 2, 3 }), isFalse());
        assertThat(GenomeDescriptorSet.CloseGenome.test(resultList, new int[] { 1, 0, 4, 2, 3 }), isFalse());
        assertThat(resultList.get(0).isSameGenome(resultList.get(4)), isFalse());
        assertThat(resultList.get(1).isSameGenome(resultList.get(2)), isTrue());
    }

}
