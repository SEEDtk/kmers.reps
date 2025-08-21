/**
 *
 */
package org.theseed.proteins.kmers.reps;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;

import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.roles.RoleUtilities;
import org.theseed.sequence.GenomeDescriptor;

/**
 * @author Bruce Parrello
 *
 */
public class TestRolePatterns {

    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(TestRolePatterns.class);


    @Test
    public void testPatterns() throws IOException {
        File gFile = new File("data", "202462.4.gto");
        Genome genome = new Genome(gFile);
        int pheCount = 0;
        int ssuCount = 0;
        for (Feature feat : genome.getFeatures()) {
            String function = RoleUtilities.commentFree(feat.getPegFunction());
            Matcher m = GenomeDescriptor.SEED_PROTEIN.matcher(function);
            if (m.matches()) {
                log.info("PHES in {}: {}", feat.getId(), function);
                pheCount++;
            } else {
                m = RoleUtilities.SSU_R_RNA.matcher(function);
                if (m.matches()) {
                    log.info("SSU in {}: {}", feat.getId(), function);
                    ssuCount++;
                }
            }
        }
        assertThat(pheCount, equalTo(1));
        assertThat(ssuCount, greaterThan(0));
    }

}
