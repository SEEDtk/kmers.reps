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
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.proteins.RoleMatrix;
import org.theseed.utils.FloatList;
import org.theseed.utils.ICommand;

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
 * -h	display the parameters
 * -v	show progress on STDERR
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
public class RolesProcessor implements ICommand {

    // FIELDS
    /** role profiles input stream */
    Scanner inStream;
    /** list of target completeness fractions */
    FloatList targets;

    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name = "-h", aliases = { "--help" }, help = true)
    protected boolean help;

    /** TRUE if we want progress messages */
    @Option(name = "-v", aliases = { "--verbose", "--debug" }, usage = "display progress on STDERR")
    protected boolean debug;

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
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        this.inFile = null;
        this.targets = new FloatList(new double[] { 0.90, 0.80 });
        this.commonOccurrence = 0.95;
        this.minRoles = 20;
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                if (this.inFile == null) {
                    this.inStream = new Scanner(System.in);
                    if (this.debug) System.err.println("Reading from standard input.");
                } else {
                    FileInputStream fileStream = new FileInputStream(this.inFile);
                    this.inStream = new Scanner(fileStream);
                    if (this.debug) System.err.println("Reading from " + this.inFile);
                }
                // Denote we can run this process.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        } catch (IOException e) {
            e.printStackTrace(System.err);
        }
        return retVal;
    }

    @Override
    public void run() {
        if (this.debug) System.err.println("Processing input.");
        // Set the delimiter to tab or new-line.
        this.inStream.useDelimiter("\t|[\r\n]+");
        // Loop until we reach end-of-file.  We process one profile group at a time.
        while (this.inStream.hasNext()) {
            // The current record has four fields: ID, name, score, and sequence.
            String groupId = this.inStream.next();
            String groupName = this.inStream.next();
            int groupScore = this.inStream.nextInt();
            String groupSeq = this.inStream.next();
            if (this.debug)
                System.err.println("Processing group " + groupId + ": " + groupName);
            // Create the role processor for this group.
            RoleMatrix roleMtx = new RoleMatrix(100, 1000);
            int gCount = 0;
            for (String g = this.inStream.next(); ! g.contentEquals("//"); g = this.inStream.next()) {
                Collection<String> roles = Arrays.asList(StringUtils.split(this.inStream.next(), ','));
                roleMtx.register(g, roles);
                gCount++;
            }
            if (this.debug) System.err.println(gCount + " genomes in group.");
            // We are searching for a good role set.  We try marker roles at each of the
            // completeness levels in the float list, and if none of them work, we fall back
            // to common roles.  The role set found is put in here.
            Collection<String> roleSet = new ArrayList<String>();
            // Start with the first target value.
            targets.reset();
            // Loop through the targets.
            while (this.targets.hasNext() && roleSet.size() < this.minRoles) {
                double target = targets.next();
                if (this.debug) System.err.println("Computing marker role set for target " + target);
                roleSet = roleMtx.getMarkerRoles(target);
            }
            if (roleSet.size() < this.minRoles) {
                // Here we failed to find a role set, so we look for common roles.
                if (this.debug) System.err.println("Computing common role set at threshold " +
                        this.commonOccurrence);
                roleSet = roleMtx.getCommonRoles(this.commonOccurrence);
            }
            if (this.debug) System.err.println(roleSet.size() + " roles found for " + groupId);
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
