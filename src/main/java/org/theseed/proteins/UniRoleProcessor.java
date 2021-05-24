/**
 *
 */
package org.theseed.proteins;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command produces a list of the singly-occurring roles for each genome in a set of input genomes.
 * Only the roles found in a specified role map are used.
 *
 * The report is written to the standard output.  Each output line consists of a role ID and the number of
 * genomes in which it occurs singly.  The output is sorted by most frequent to least frequent.
 *
 * The positional parameters are the name of the role map file and the name of the input file or directory.  The input
 * is a genome source-- either a master directory, a GTO directory (the default), or a file of PATRIC genome IDs.  The
 * role map is a tab-delimited file with role IDs in the first column and role descriptions in the third column.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output for report (if not STDOUT)
 * -m	minimum percent of genomes that must contain the role for it to be included in the output; default 80
 *
 * --source		specify the type of input (PATRIC, DIR, or MASTER); default DIR
 *
 * @author Bruce Parrello
 */
public class UniRoleProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UniRoleProcessor.class);
    /** role definitions */
    private RoleMap roleMap;
    /** genome source */
    private GenomeSource genomes;
    /** output stream */
    private OutputStream output;
    /** count of role occurrences */
    private CountMap<String> roleCounts;
    /** map of role IDs to length statistics */
    private Map<String, SummaryStatistics> roleStats;
    /** map of role IDs to protein families */
    private Map<String, Set<String>> familyMap;
    /** output threshold */
    private int minCount;
    /** estimated number of roles for map creation */
    private static final int ESTIMATED_ROLES = 5000;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome input (master genome directory, GTO directory, patric ID file)")
    private GenomeSource.Type inType;

    /** output file (if not STDOUT) */
    @Option(name = "--output", aliases = { "-o" }, usage = "output file (if not STDOUT)")
    private File outFile;

    /** minimum percent of genomes that must contain a role */
    @Option(name = "--min", aliases = { "-m" }, usage = "minimum percent of genomes that must contain a role")
    private double minPercent;

    /** role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "role definition file", required = true)
    private File roleFile;

    /** input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome file or directory", required = true)
    private File inDir;


    @Override
    protected void setDefaults() {
        this.inType = GenomeSource.Type.DIR;
        this.outFile = null;
        this.minPercent = 80.0;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the role file.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " not found or unreadable.");
        this.roleMap = RoleMap.load(this.roleFile);
        log.info("{} role definitions loaded from {}.", this.roleFile);
        // Set up the output.
        if (this.outFile == null) {
            log.info("Report will be written to standard output.");
            this.output = System.out;
        } else {
            log.info("Report will be written to {}.", this.outFile);
            this.output = new FileOutputStream(this.outFile);
        }
        // Verify the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " not found.");
        this.genomes = this.inType.create(this.inDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.inDir);
        // Compute the output threshold.
        if (this.minPercent > 100.0)
            throw new ParseFailureException("Minimum percent must be <= 100");
        this.minCount = (int) (this.minPercent * this.genomes.size() + 99) / 100;
        log.info("Output threshold is {}% ({} occurrences).", this.minPercent, this.minCount);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // Create the role maps.
            this.roleCounts = new CountMap<String>();
            this.roleStats = new HashMap<String, SummaryStatistics>(ESTIMATED_ROLES);
            this.familyMap = new HashMap<String, Set<String>>(ESTIMATED_ROLES);
            // Loop through the genomes.
            int count = 0;
            for (Genome genome : genomes) {
                count++;
                log.info("Processing genome {} of {}: {}.", count, genomes.size(), genome);
                // Count the role occurrences.
                CountMap<String> gRoleCounts = new CountMap<String>();
                CountMap<String> gRoleLengths = new CountMap<String>();
                for (Feature feat : genome.getPegs()) {
                    for (Role role : feat.getUsefulRoles(this.roleMap)) {
                        String roleId = role.getId();
                        gRoleCounts.count(roleId);
                        gRoleLengths.setCount(roleId, feat.getProteinLength());
                        this.updateFamilies(roleId, feat);
                    }
                }
                // Now output the roles with only one feature.  For each role we must also add its length
                // to the statistics for that role.
                int rCount = 0;
                for (CountMap<String>.Count counter : gRoleCounts.counts()) {
                    if (counter.getCount() == 1) {
                        String roleId = counter.getKey();
                        this.roleCounts.count(roleId);
                        SummaryStatistics stats = this.roleStats.computeIfAbsent(roleId, x -> new SummaryStatistics());
                        stats.addValue(gRoleLengths.getCount(roleId));
                        rCount++;
                    }
                }
                log.info("{} singleton roles found for {}.", rCount, genome);
            }
            if (this.roleCounts.size() == 0)
                log.error("No singleton roles found.");
            else {
                log.info("Producing output. {} distinct singleton roles found.", this.roleCounts.size());
                try (PrintWriter writer = new PrintWriter(this.output)) {
                    writer.println("role\tcount\tname\tmean_len\tsdev_len\tfamilies");
                    // Loop through the roles of interest.
                    for (CountMap<String>.Count counter : this.roleCounts.sortedCounts()) {
                        int rCount = counter.getCount();
                        if (rCount >= this.minCount) {
                            String roleId = counter.getKey();
                            // Get the protein families.
                            String families = StringUtils.join(this.familyMap.get(roleId), ", ");
                            // We need to compute the mean and standard deviation.
                            SummaryStatistics stats = this.roleStats.get(roleId);
                            writer.format("%s\t%d\t%s\t%8.2f\t%8.2f\t%s%n", roleId, rCount, this.roleMap.getName(roleId),
                                    stats.getMean(), stats.getStandardDeviation(), families);
                        }
                    }
                }
            }
        } finally {
            // Insure we close the output stream.  This really only matters if we are
            // being called as a service by a client program, since in that case we
            // want the file closed immediately, rather than waiting for the process to
            // finish.
            if (this.outFile != null)
                this.output.close();
        }
    }

    /**
     * Add the global pattyfam and figfam for a role to the family map.
     *
     * @param roleId	ID of the role of interest
     * @param feat		feature containing the role
     */
    private void updateFamilies(String roleId, Feature feat) {
        Set<String> families = this.familyMap.computeIfAbsent(roleId, x -> new TreeSet<String>());
        if (feat.getFigfam() != null) families.add(feat.getFigfam());
        if (feat.getPgfam() != null) families.add(feat.getPgfam());
    }

}
