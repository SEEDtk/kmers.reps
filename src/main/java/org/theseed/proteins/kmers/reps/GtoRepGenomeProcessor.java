/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;

/**
 * This command will process all the GTOs in a directory and determine the closest representative genome to each.
 * The positional parameters are the name of the representative-genome database and the name of the genome
 * directory.  The report will be on the standard output.  If a genome does not have a seed protein it will be
 * omitted from the report.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more progress messages on the log
 *
 * @author Bruce Parrello
 *
 */
public class GtoRepGenomeProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoRepGenomeProcessor.class);
    /** representative-genome database */
    private RepGenomeDb repDb;
    /** role map for seed protein */
    private RoleMap seedMap;

    // COMMAND-LINE OPTIONS

    /** representative-genome database */
    @Argument(index = 0, metaVar = "repDb.ser", usage = "name of representaitve-genome database file")
    private File repDbFile;

    /** input genome directory */
    @Argument(index = 1, metaVar = "gtoDir", usage = "input GTO directory")
    private File gtoDir;


    @Override
    protected void setDefaults() {
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Representative-genome database " + this.repDbFile + " not found or unreadable.");
        if (! this.gtoDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.gtoDir + " not found or invalid.");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Load the rep-genome database.
        log.info("Loading representative-genome database from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
        // Start the output report.
        System.out.println("genome_id\tgenome_name\trep_id\trep_name\tsimilarity\tdistance");
        // Create a role map from the seed protein.
        this.seedMap = new RoleMap();
        this.seedMap.findOrInsert(this.repDb.getProtName());
        // Now read through the genome directory.
        log.info("Scanning genome directory {}.", this.gtoDir);
        GenomeDirectory genomes = new GenomeDirectory(this.gtoDir);
        log.info("{} genomes found.", genomes.size());
        for (Genome genome : genomes) {
            log.info("Processing genome {}.", genome);
            // Search for the seed protein.  We keep the longest.
            String seedFound = "";
            for (Feature feat : genome.getPegs()) {
                String protein = feat.getProteinTranslation();
                // If this protein is longer, parse the functional assignment.
                if (protein.length() > seedFound.length()) {
                    List<Role> seeds = feat.getUsefulRoles(this.seedMap);
                    if (seeds.size() > 0)
                        seedFound = protein;
                }
            }
            // If we have a seed, find its closest representative.
            if (seedFound.isEmpty())
                log.warn("{} does not have a seed protein.", genome);
            else {
                RepGenomeDb.Representation result = this.repDb.findClosest(seedFound);
                String repId;
                String repName;
                if (result.getSimilarity() == 0) {
                    // Here there is nothing close.
                    repId = "";
                    repName = "";
                } else {
                    // Compute the ID and name of the closest rep.
                    repId = result.getGenomeId();
                    repName = result.getRepresentative().getName();
                }
                // Write the output line.
                System.out.format("%s\t%s\t%s\t%s\t%d\t%8.4f%n", genome.getId(), genome.getName(), repId, repName,
                        result.getSimilarity(), result.getDistance());
            }
        }
    }
}
