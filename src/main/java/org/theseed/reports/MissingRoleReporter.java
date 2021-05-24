/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;

import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleType;

/**
 * This object produces reports for the missing-role finder.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MissingRoleReporter extends BaseWritingReporter {

    /**
     * This enum determines the type of report output.
     */
    public static enum Type {
        FULL {
            @Override
            public MissingRoleReporter create(PrintWriter writer) {
                return new TaxonBasedMissingRoleReporter(writer, RoleType.values());
            }
        }, MISSING {
            @Override
            public MissingRoleReporter create(PrintWriter writer) {
                return new TaxonBasedMissingRoleReporter(writer, RoleType.MISSING);
            }
        }, OVER_PRESENT {
            @Override
            public MissingRoleReporter create(PrintWriter writer) {
                return new TaxonBasedMissingRoleReporter(writer, RoleType.OVER_PRESENT);
            }
        }, BAD {
            @Override
            public MissingRoleReporter create(PrintWriter writer) {
                return new TaxonBasedMissingRoleReporter(writer, RoleType.MISSING, RoleType.OVER_PRESENT);
            }
        };

        /**
         * @return a missing-role reporter of the specified type
         *
         * @param writer	print writer to receive the output
         */
        public abstract MissingRoleReporter create(PrintWriter writer);
    }

    /**
     * Construct a new missing-role reporter.
     *
     * @param writer	print writer to receive the output
     */
    public MissingRoleReporter(PrintWriter writer) {
        super(writer);
    }

    /**
     * Initialize the report.
     */
    public abstract void openReport();

    /**
     * Record the number of occurrences of a role in a genome.  Every role will be reported for every genome.
     *
     * @param role		descriptor for the role of interest
     * @param count		number of times it occurred in the genome
     * @param genome	genome of interest
     */
    public abstract void recordRole(Role role, int count, Genome genome);

    /**
     * Finish the report.  If the output is sorted, it will all be produced here.
     */
    public abstract void finish();

}
