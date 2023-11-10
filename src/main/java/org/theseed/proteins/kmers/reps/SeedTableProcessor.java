/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseReportProcessor;

/**
 * This command creates a 4-column table from an input genome source.  The four-column table
 * is based on a specific seed protein, and contains for each genome (1) the genome ID, (2) the
 * genome name, (3) the seed-protein amino acid sequence, and (4) the seed-protein DNA sequence.
 *
 * The positional parameter is the name of the input genome source and the ID of the seed protein.
 * All the genomes in the source will be processed.  If the identified seed protein is not present
 * a warning will be issued and the genome will be omitted in the output file.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 *
 * --source		type of genome source (default DIR)
 * --roles		role definition table (default "roles.in.subsystems" in the current directory)
 *
 * @author Bruce Parrello
 *
 */
public class SeedTableProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SeedTableProcessor.class);
    /** input genome source */
    private GenomeSource source;
    /** role definition table (contains only the seed role) */
    private RoleMap roleMap;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "genome source type")
    private GenomeSource.Type sourceType;

    /** role definition file */
    @Option(name = "--roles", metaVar = "roles.in.subsystems", usage = "role definition file")
    private File roleFile;

    /** genome source directory */
    @Argument(index = 0, metaVar = "genomeDir", usage = "genome source directory", required = true)
    private File genomeDir;

    /** seed protein role ID */
    @Argument(index = 1, metaVar = "roleId", usage = "seed protein role ID", required = true)
    private String roleId;

    @Override
    protected void setReporterDefaults() {
        this.roleFile = new File(System.getProperty("user.dir"), "roles.in.subsystems");
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Insure the genome source exists.
        if (! this.genomeDir.exists())
            throw new FileNotFoundException("Genome source " + this.genomeDir + " not found.");
        // Get the role definition file and extract the identified seed role into the role map.
        this.roleMap = new RoleMap();
        try (TabbedLineReader roleStream = new TabbedLineReader(this.roleFile, 3)) {
            for (TabbedLineReader.Line line : roleStream) {
                String inRole = line.get(0);
                if (inRole.contentEquals(this.roleId))
                    this.roleMap.addRole(inRole, line.get(2));
            }
            // Insure we found the role.
            if (this.roleMap.size() < 1)
                throw new ParseFailureException("Role \"" + this.roleId + "\" not found in role file "
                        + this.roleFile.toString() + ".");
        }
        // Now open the genome source for input.
        this.source = this.sourceType.create(this.genomeDir);
        log.info("{} genomes found in source {} of type {}.", this.source.size(), this.genomeDir,
                this.sourceType);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the output headings.
        writer.format("genome_id\tgenome_name\t%s_protein\t%s_dna%n",
                this.roleId, this.roleId);
        // This will track our statistics.
        int missing = 0;
        int found = 0;
        // Loop through the genomes.
        for (Genome genome : this.source) {
            log.info("Processing {}.", genome);
            // Loop through the features.  We keep the longest seed protein.
            String bestProtein = "";
            String bestDna = "";
            for (Feature feat : genome.getPegs()) {
                if (feat.isInteresting(this.roleMap)) {
                    // Here we have found a feature.
                    String thisProtein = feat.getProteinTranslation();
                    if (thisProtein.length() > bestProtein.length()) {
                        // This protein is better than the one we have, so we use it.
                        bestDna = genome.getDna(feat.getLocation());
                        bestProtein = thisProtein;
                    }
                }
            }
            // If we did not find a good protein, it's a warning.
            if (bestProtein.isEmpty()) {
                log.warn("No {} protein found for {}.", this.roleId, genome);
                missing++;
            } else {
                // Otherwise, write a record.
                writer.format("%s\t%s\t%s\t%s%n", genome.getId(), genome.getName(),
                        bestProtein, bestDna);
                found++;
            }
        }
        // All done.  Write the stats.
        log.info("{} proteins found, {} missing.", found, missing);
    }

}
