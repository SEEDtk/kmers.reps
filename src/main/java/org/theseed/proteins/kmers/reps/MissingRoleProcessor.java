/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.reports.MissingRoleReporter;
import org.theseed.utils.BaseReportProcessor;

/**
 * This command compares a list of universal roles (produced by UniRoleProcessor) to a directory of genomes.
 * It will output a table indicating whether each role is present, absent, or over-present in each genome.  The
 * table will be sorted by role ID followed by taxonomy.
 *
 * The positional parameters are the name of the role definition file (roles.in.subsystems), the name of the
 * universal role file, and finally the name of the genome input source.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	name of the report output file (if not STDOUT)
 * -t	type of input source (PATRIC, DIR, MASTER)
 *
 * --format output report format
 *
 * @author Bruce Parrello
 *
 */
public class MissingRoleProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MissingRoleProcessor.class);
    /** map of universal roles */
    private RoleMap roleMap;
    /** genome input source */
    private GenomeSource source;
    /** reporting object */
    private MissingRoleReporter reporter;

    // COMMAND-LINE OPTIONS

    /** input source type */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** output report format */
    @Option(name = "--format", usage = "output report format")
    private MissingRoleReporter.Type outputType;

    /** role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "current role definition file", required = true)
    private File roleFile;

    /** universal role input file */
    @Argument(index = 1, metaVar = "sours.tbl", usage = "input file containing universal role IDs and names", required = true)
    private File sourFile;

    /** input genome source */
    @Argument(index = 2, metaVar = "genomeDir", usage = "input genome source (file or directory)")
    private File genomeDir;

    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.outputType = MissingRoleReporter.Type.FULL;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify that the input role files exist.
        if (! this.sourFile.canRead())
            throw new FileNotFoundException("Universal role file " + this.sourFile + " is not found or unreadable.");
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or unreadable.");
        // Connect to the genome source.
        this.source = this.sourceType.create(genomeDir);
        log.info("{} genomes found in {} source {}.", this.source.size(), this.sourceType, this.genomeDir);
        // We need to create the role map.  This is a two-stage process. First, we read the universal roles.  Next,
        // we create the real role map by keeping only the role definitions for the universals.  The key thing here
        // is that some role IDs have aliases, and will be added multiple times.
        Set<String> uniRoles = new HashSet<String>(300);
        try (TabbedLineReader reader = new TabbedLineReader(this.sourFile)) {
            for (TabbedLineReader.Line line : reader)
                uniRoles.add(line.get(0));
        }
        // Note that the role definition file has no header.
        this.roleMap = new RoleMap();
        try (TabbedLineReader reader = new TabbedLineReader(this.roleFile, 3)) {
            for (TabbedLineReader.Line line : reader) {
                String roleId = line.get(0);
                if (uniRoles.contains(roleId)) {
                    this.roleMap.addRole(roleId, line.get(2));
                }
            }
        }
        log.info("{} universal roles found in {} with {} aliases included.", this.roleMap.size(), this.sourFile,
                this.roleMap.fullSize() - this.roleMap.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Start the report.
        this.reporter = this.outputType.create(writer);
        this.reporter.openReport();
        // Loop through the genomes.  We count the number of genomes processed and the number of
        // missing roles.  We also compute statistics for the number of roles per genome.
        int count = 0;
        DescriptiveStatistics stats = new DescriptiveStatistics();
        int missingCount = 0;
        for (Genome genome : this.source) {
            count++;
            log.info("Processing genome {} of {}:  {}", count, this.source.size(), genome);
            // Loop through the features, collecting universal role counts.
            CountMap<String> roleCounts = new CountMap<String>();
            for (Feature feat : genome.getPegs()) {
                // Note that the map contains only universal roles, so ONLY universals will
                // be found here.
                for (Role role : feat.getUsefulRoles(this.roleMap))
                    roleCounts.count(role.getId());
            }
            log.info("{} distinct universal roles found in {}.", roleCounts.size(), genome);
            stats.addValue((double) roleCounts.size());
            // Output the roles for this genome.
            for (Role role : this.roleMap) {
                int rCount = roleCounts.getCount(role.getId());
                if (rCount == 0) missingCount++;
                this.reporter.recordRole(role, rCount, genome);
            }
        }
        // Finish the report.
        this.reporter.finish();
        // Log the statistics.
        log.info("{} roles were missing.", missingCount);
        log.info("Roles/genome minimum = {}, maximum = {}, mean = {}, median = {}.", stats.getMin(), stats.getMax(),
                stats.getMean(), stats.getPercentile(50.0));
    }

}
