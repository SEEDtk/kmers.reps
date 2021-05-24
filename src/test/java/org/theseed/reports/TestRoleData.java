/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import org.theseed.genome.Genome;
import org.theseed.proteins.Role;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * @author Bruce Parrello
 *
 */
public class TestRoleData {

    private static final String[] GTOS = new String[] { "103621.4.gto", "1036677.3.gto", "1036778.3.gto", "1280.21158.gto",
            "1280.25101.gto", "1313.16015.gto", "202462.4.gto" };

    @Test
    public void testRoleData() throws IOException {
        Genome[] genomes = new Genome[GTOS.length];
        for (int i = 0; i < GTOS.length; i++) {
            File gFile = new File("data", GTOS[i]);
            genomes[i] = new Genome(gFile);
        }
        List<RoleData> roles = new ArrayList<RoleData>(GTOS.length * 2);
        Role pheS = new Role("PhenTrnaSyntAlph", "Phenylalanyl-tRNA synthetase alpha chain");
        Role lrna = new Role("LsuRiboProt", "LSU ribosomal protein");
        for (int i = 0; i < GTOS.length; i++) {
            roles.add(new RoleData(pheS, i % 3, genomes[i]));
            roles.add(new RoleData(lrna, 2 - (i % 3), genomes[i]));
        }
        Collections.sort(roles);
        assertThat(roles.get(0).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(0).getGenomeId(), equalTo("1036778.3"));
        assertThat(roles.get(1).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(1).getGenomeId(), equalTo("1280.21158"));
        assertThat(roles.get(2).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(2).getGenomeId(), equalTo("1280.25101"));
        assertThat(roles.get(3).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(3).getGenomeId(), equalTo("1313.16015"));
        assertThat(roles.get(4).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(4).getGenomeId(), equalTo("103621.4"));
        assertThat(roles.get(5).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(5).getGenomeId(), equalTo("202462.4"));
        assertThat(roles.get(6).getRoleId(), equalTo("LsuRiboProt"));
        assertThat(roles.get(6).getGenomeId(), equalTo("1036677.3"));
        assertThat(roles.get(7).getRoleId(), equalTo("PhenTrnaSyntAlph"));
        assertThat(roles.get(7).getGenomeId(), equalTo("1036778.3"));
    }

}
