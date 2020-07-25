/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Scanner;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Option;
import org.theseed.proteins.RoleMatrix;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.FloatList;

/**
 * This class processes a universal-role file to create a rep-genome based completeness engine.
 * For each group, we need the group ID, the similarity score (0 for taxonomy-based groups), the
 * seed protein (empty for taxonomy-based groups), and a comma-delimited list of the universal roles.
 * Together with a file mapping role names to IDs, this is enough to build a completeness processor.
 *
 * The input file consists of one or more groups.  Each group has a header, one record per genome,
 * and a trailer.  Every record is tab-delimited.  The header contains the group ID, group name,
 * similarity score, and seed protein sequence.  There will be one group for leftover genomes with
 * a name of "root".  Each data record contains a genome ID and a comma-delimited role list.
 * The trailer is simply the string "//".
 *
 * The command-line options are as follows.
 *
 * -i	specifies a file name to use for the standard input
 * -m	minimum number of acceptable roles for a completeness set
 *
 * --target		target completeness fractions for computing marker roles, coded as a comma-delimited
 * 				list; if a particular fraction doesn't yield enough roles, the next one will be used;
 * 				if none work, then common roles will be sought
 * --common		minimum commonality fraction for computing common roles
 *
 * @author Bruce Parrello
 *
 */
public class RolesProcessor extends BaseProcessor {

    // FIELDS
    /** role profiles input stream */
    Scanner inStream;
    /** list of target completeness fractions */
    FloatList targets;

    // COMMAND-LINE OPTIONS

    /** input file (if not using STDIN) */
    @Option(name = "--input", aliases = { "-i" }, usage = "input file (if not STDIN)")
    private File inFile;

    /** target completeness fractions */
    @Option(name = "--target", usage="comma-delimited list of target completeness fractions for computing marker roles")
    private void setTargets(String targetString) {
        this.targets = new FloatList(targetString);
    }

    /** minimum commonality fraction */
    @Option(name="--common", usage="target commonality fraction for computing common roles")
    private double commonOccurrence;

    /** minimum number of roles for an acceptable role set */
    @Option(name="-m", aliases = { "--minRoles", "--min" }, usage="minimum number of roles for a set")
    private int minRoles;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.targets = new FloatList(new double[] { 0.90, 0.80 });
        this.commonOccurrence = 0.95;
        this.minRoles = 20;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (this.inFile == null) {
            this.inStream = new Scanner(System.in);
            log.info("Reading from standard input.");
        } else {
            FileInputStream fileStream = new FileInputStream(this.inFile);
            this.inStream = new Scanner(fileStream);
            log.info("Reading from {}.", this.inFile);
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        log.info("Processing input.");
        // Set the delimiter to tab or new-line.
        this.inStream.useDelimiter("\t|[\r\n]+");
        // Loop until we reach end-of-file.  We process one profile group at a time.
        while (this.inStream.hasNext()) {
            // The current record has four fields: ID, name, score, and sequence.
            String groupId = this.inStream.next();
            String groupName = this.inStream.next();
            int groupScore = this.inStream.nextInt();
            String groupSeq = this.inStream.next();
            log.info("Processing group {}: {}.", groupId, groupName);
            // Create the role processor for this group.
            RoleMatrix roleMtx = new RoleMatrix(100, 1000);
            int gCount = 0;
            for (String g = this.inStream.next(); ! g.contentEquals("//"); g = this.inStream.next()) {
                Collection<String> roles = Arrays.asList(StringUtils.split(this.inStream.next(), ','));
                roleMtx.register(g, roles);
                gCount++;
            }
            log.info("{} genomes in group.", gCount);
            // We are searching for a good role set.  We try marker roles at each of the
            // completeness levels in the float list, and if none of them work, we fall back
            // to common roles.  The role set found is put in here.
            Collection<String> roleSet = new ArrayList<String>();
            // Start with the first target value.
            targets.reset();
            // Loop through the targets.
            double target = 0.0;
            while (this.targets.hasNext() && roleSet.size() < this.minRoles) {
                target = targets.next();
                log.debug("Computing marker role set for target {}.");
                roleSet = roleMtx.getMarkerRoles(target);
            }
            if (roleSet.size() < this.minRoles) {
                // Here we failed to find a role set, so we look for common roles.
                log.info("Computing common role set at threshold {}.", this.commonOccurrence);
                roleSet = roleMtx.getCommonRoles(this.commonOccurrence);
                // Count the genomes that fail the lowest target.
                int fCount = 0;
                for (String genome : roleMtx.genomes()) {
                    if (roleMtx.completeness(roleSet, genome) < target) fCount++;
                }
                log.info("{} genomes failed the lowest target.", fCount);
            }
            log.info("{} roles found for {}.", roleSet.size(), groupId);
            // Write out the group header.
            System.out.format("%s\t%d\t%s\t%s%n", groupId, groupScore, groupName, groupSeq);
            // Write out the roles.
            for (String role : roleSet) {
                System.out.println("   " + role);
            }
            // Write out the trailer.
            System.out.println("//");
        }
    }

}
