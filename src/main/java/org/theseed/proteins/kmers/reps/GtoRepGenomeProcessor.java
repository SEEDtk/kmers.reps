/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command will process all the genomes from a source and determine the closest representative genome to each.
 * The positional parameters are the name of the representative-genome database and the name of the genome
 * source (this can be a master directory, a file of PATRIC genome IDs, or a normal GTO directory).  The report will be
 * on the standard output.  If a genome does not have a seed protein it will be omitted from the report.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more progress messages on the log
 * -t	type of genome directory source (MASTER, DIR, PATRIC)
 *
 * --filter		if specified, a tab-delimited file of genome IDs; only the genome IDs listed in the first column
 * 				will be processed)
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
    /** input genome source */
    private GenomeSource source;
    /** genome ID filter */
    private Set<String> idFilter;

    // COMMAND-LINE OPTIONS

    /** genome source type */
    @Option(name = "-t", aliases = { "--type" }, usage = "genome source type")
    private GenomeSource.Type type;

    /** filter file */
    @Option(name = "--filter", metaVar = "idFile.tbl", usage = "if specified, a tab-delimited file of genome IDs to select")
    private File filterFile;

    /** representative-genome database */
    @Argument(index = 0, metaVar = "repDb.ser", usage = "name of representaitve-genome database file")
    private File repDbFile;

    /** input genome directory */
    @Argument(index = 1, metaVar = "gtoDir", usage = "input GTO source")
    private File gtoDir;


    @Override
    protected void setDefaults() {
        this.type = GenomeSource.Type.DIR;
        this.filterFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Representative-genome database " + this.repDbFile + " not found or unreadable.");
        if (! this.gtoDir.exists())
            throw new FileNotFoundException("Input source " + this.gtoDir + " not found.");
        // If there is a filter file, create the filter set.
        if (this.filterFile == null) {
            log.info("All input genomes will be processed.");
            this.idFilter = null;
        } else {
            this.idFilter = TabbedLineReader.readSet(this.filterFile, "1");
            log.info("{} genomes specified in filter file {}.", this.idFilter.size(), this.filterFile);
        }
        // Load the genome source.
        this.source = this.type.create(this.gtoDir);
        log.info("{} genomes in source {}.", this.source.size(), this.gtoDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Load the rep-genome database.
        log.info("Loading representative-genome database from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
        // Start the output report.
        System.out.println("genome_id\tgenome_name\trep_id\trep_name\tsimilarity\tdistance\toutlier");
        // Create a role map from the seed protein.
        this.seedMap = new RoleMap();
        this.seedMap.findOrInsert(this.repDb.getProtName());
        // Now read through the genome directory.
        log.info("Scanning genome directory {}.", this.gtoDir);
        int filterCount = 0;
        int skipCount = 0;
        int outCount = 0;
        for (String genomeId : this.source.getIDs()) {
            if (this.idFilter != null && ! this.idFilter.contains(genomeId))
                filterCount++;
            else {
                Genome genome = this.source.getGenome(genomeId);
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
                if (seedFound.isEmpty()) {
                    log.warn("{} does not have a seed protein.", genome);
                    skipCount++;
                } else {
                    RepGenomeDb.Representation result = this.repDb.findClosest(seedFound);
                    String repId;
                    String repName;
                    String outlierFlag = "";
                    if (result.getSimilarity() == 0) {
                        // Here there is nothing close.
                        repId = "";
                        repName = "";
                        outCount++;
                        outlierFlag = "*";
                    } else {
                        // Compute the ID and name of the closest rep.
                        repId = result.getGenomeId();
                        repName = result.getRepresentative().getName();
                        if (result.getSimilarity() < repDb.getThreshold()) {
                            outCount++;
                            outlierFlag = "*";
                        }
                    }
                    // Write the output line.
                    System.out.format("%s\t%s\t%s\t%s\t%d\t%8.4f\t%s%n", genome.getId(), genome.getName(), repId, repName,
                            result.getSimilarity(), result.getDistance(), outlierFlag);
                }
            }
        }
        log.info("{} outliers found, {} skipped due to filtering, {} skipped due to missing seed protein.",
                outCount, filterCount, skipCount);
    }
}
