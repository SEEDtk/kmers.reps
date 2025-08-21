/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleType;

/**
 * This is an exhaustive missing-role report that allows present, absent, and/or over-present  roles.  The report is sorted by role ID
 * and then taxonomy so that the clade distribution of the missing roles can be seen.
 *
 * @author Bruce Parrello
 *
 */
public class TaxonBasedMissingRoleReporter extends MissingRoleReporter {

    // FIELDS
    /** output role data queue */
    private SortedSet<RoleData> dataLines;
    /** type of roles to report */
    private final Set<RoleType> allowed;

    /**
     * Create a report for the specified role types.
     *
     * @param writer	output print writer
     * @param roleTypes	acceptable role types
     */
    public TaxonBasedMissingRoleReporter(PrintWriter writer, RoleType... roleTypes) {
        super(writer);
        this.allowed = new TreeSet<>(Arrays.asList(roleTypes));
    }

    @Override
    public void openReport() {
        // Write the header.
        this.println(RoleData.getHeader());
        // Create the data line tree.
        this.dataLines = new TreeSet<>();
    }

    @Override
    public void recordRole(Role role, int count, Genome genome) {
        // Keep this data line if it's a type we want.
        RoleType roleType = RoleType.compute(count);
        if (this.allowed.contains(roleType))
            this.dataLines.add(new RoleData(role, count, genome));
    }

    @Override
    public void finish() {
        // Unspool the data lines.
        for (RoleData dataLine : this.dataLines)
            this.println(dataLine.output());
    }

}
