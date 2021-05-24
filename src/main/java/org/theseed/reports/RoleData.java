/**
 *
 */
package org.theseed.reports;

import org.theseed.genome.Genome;
import org.theseed.proteins.Role;

/**
 * This object contains a single line of data for a missing-roles report. It is sorted by role ID and then taxonomy.
 */
public class RoleData implements Comparable<RoleData> {

    // FIELDS
    /** ID of the role */
    private String roleId;
    /** descriptive name of the role */
    private String roleName;
    /** number of occurrences */
    private int count;
    /** ID of the genome */
    private String genomeId;
    /** taxonomy string of the genome */
    private String taxString;
    /** lineage array, used for sorting */
    private int[] lineage;

    /**
     * Create a role data descriptor.
     *
     * @param role		role being counted
     * @param count		number of role occurrences
     * @param genome	genome containing/missing the role
     */
    public RoleData(Role role, int count, Genome genome) {
        this.roleId = role.getId();
        this.roleName = role.getName();
        this.count = count;
        this.genomeId = genome.getId();
        this.taxString = genome.getTaxString();
        this.lineage = genome.getLineage();
    }

    @Override
    public int compareTo(RoleData o) {
        int retVal = roleId.compareTo(o.roleId);
        for (int i = 0; retVal == 0 && i < this.lineage.length; i++) {
            if (i >= o.lineage.length)
                retVal = 1;
            else
                retVal = this.lineage[i] - o.lineage[i];
        }
        if (retVal == 0) {
            if (this.lineage.length < o.lineage.length)
                retVal = -1;
            else
                retVal = this.genomeId.compareTo(o.genomeId);
        }
        return retVal;
    }

    /**
     * @return the header for the report
     */
    public static String getHeader() {
        return "role_id\trole_name\tcount\tgenome_id\ttaxonomy";
    }

    /**
     * @return the output for this data line
     */
    public String output() {
        return String.format("%s\t%s\t%d\t%s\t%s", this.roleId, this.roleName, this.count, this.genomeId, this.taxString);
    }

    /**
     * @return the role ID
     */
    public String getRoleId() {
        return this.roleId;
    }

    /**
     * @return the genome ID
     */
    public String getGenomeId() {
        return this.genomeId;
    }

}
