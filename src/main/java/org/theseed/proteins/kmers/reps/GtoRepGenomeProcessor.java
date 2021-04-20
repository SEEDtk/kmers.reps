/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.genome.iterator.GenomeTargetType;
import org.theseed.genome.iterator.IGenomeTarget;
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
 * -o	output file (if not STDOUT)
 *
 * --filter		if specified, a tab-delimited file of genome IDs; only the genome IDs listed in the first column
 * 				will be processed, and they will be processed in the order presented
 * --rna		if specified, genomes without a valid SSU rRNA will be skipped
 * --update		if specified, a file name; the repGenomeDB will be updated and stored in the named file
 * --save		if specified, a genome target to which outliers should be saved
 * --saveType	type of genome target (MASTER or DIR)
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
    /** genome ID filter list */
    private List<String> idList;
    /** output genome target */
    private IGenomeTarget target;
    /** output stream */
    private OutputStream outStream;


    // COMMAND-LINE OPTIONS

    /** genome source type */
    @Option(name = "-t", aliases = { "--type" }, usage = "genome source type")
    private GenomeSource.Type type;

    /** filter file */
    @Option(name = "--filter", metaVar = "idFile.tbl", usage = "if specified, a tab-delimited file of genome IDs to select")
    private File filterFile;

    /** update file for repDB */
    @Option(name = "--update", metaVar = "newRepDb.ser", usage = "optional file to contain updated repGenomeDb")
    private File updateFile;

    /** optional output genome directory */
    @Option(name = "--save", metaVar = "genomeDir", usage = "optional output directory for outlier genomes")
    private File saveDir;

    /** type of output genome directory */
    @Option(name = "--saveType", usage = "type of output directory")
    private GenomeTargetType saveType;

    /** output file (if not STDOUT) */
    @Option(name = "-o", aliases = { "--output" }, usage = "output file for report (if not STDOUT)")
    private File outFile;

    /** RNA validation flag */
    @Option(name = "--rna", usage = "if specified, genomes without a valid SSU rRNA will be skipped")
    private boolean rnaFlag;

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
        this.updateFile = null;
        this.saveType = GenomeTargetType.DIR;
        this.saveDir = null;
        this.outFile = null;
        this.rnaFlag = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Representative-genome database " + this.repDbFile + " not found or unreadable.");
        if (! this.gtoDir.exists())
            throw new FileNotFoundException("Input source " + this.gtoDir + " not found.");
        // Load the genome source.
        this.source = this.type.create(this.gtoDir);
        log.info("{} genomes in source {}.", this.source.size(), this.gtoDir);
        // If there is a filter file, create the filter set.
        if (this.filterFile == null) {
            log.info("All input genomes will be processed.");
            this.idList = new ArrayList<String>(this.source.getIDs());
        } else {
            this.idList = TabbedLineReader.readColumn(this.filterFile, "1");
            log.info("{} genomes specified in filter file {}.", this.idList.size(), this.filterFile);
        }
        // Set up the outlier output directory.
        if (this.saveDir == null)
            this.target = null;
        else {
            log.info("Outlier genomes will be output to {}.", this.saveDir);
            // Note that if the output directory does not exist we give instructions to create it.
            this.target = this.saveType.create(this.saveDir, ! this.saveDir.isDirectory());
        }
        if (this.outFile == null) {
            log.info("Output will be to the standard output.");
            this.outStream = System.out;
        } else {
            log.info("Output will be to {}.", this.outFile);
            this.outStream = new FileOutputStream(this.outFile);
        }
        // Report on the update status.
        if (this.updateFile != null)
            log.info("Outliers will be added to the database and the updated database written to {}.", this.updateFile);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        try (PrintWriter writer = new PrintWriter(this.outStream)) {
            // Load the rep-genome database.
            log.info("Loading representative-genome database from {}.", this.repDbFile);
            this.repDb = RepGenomeDb.load(this.repDbFile);
            // Start the output report.
            writer.println("genome_id\tgenome_name\trep_id\trep_name\tsimilarity\tdistance\toutlier");
            // Create a role map from the seed protein.
            this.seedMap = new RoleMap();
            this.seedMap.findOrInsert(this.repDb.getProtName());
            // Now read through the genome directory.
            log.info("Scanning genome directory {}.", this.gtoDir);
            int skipCount = 0;
            int outCount = 0;
            for (String genomeId : this.idList) {
                Genome genome = this.source.getGenome(genomeId);
                // It is completely acceptable for the genome to not be in the source.
                // We only process if it's found.
                if (genome != null) {
                    log.info("Processing genome {}.", genome);
                    // Search for the seed protein.  We keep the longest.
                    String seedFound = "";
                    String seedFid = "";
                    for (Feature feat : genome.getPegs()) {
                        String protein = feat.getProteinTranslation();
                        // If this protein is longer, parse the functional assignment.
                        if (protein.length() > seedFound.length()) {
                            List<Role> seeds = feat.getUsefulRoles(this.seedMap);
                            if (seeds.size() > 0) {
                                seedFound = protein;
                                seedFid = feat.getId();
                            }
                        }
                    }
                    // Only proceed if we have a seed protein.  Here we also check for the SSU rRNA.
                    boolean skip = false;
                    if (seedFound.isEmpty()) {
                        log.warn("{} does not have a seed protein.", genome);
                        skip = true;
                    } else if (this.rnaFlag) {
                        String ssuRRNA = genome.getSsuRRna();
                        if (ssuRRNA.isEmpty()) {
                            log.warn("{} does not have an SSU rRNA.", genome);
                            skip = true;
                        }
                    }
                    if (skip)
                        skipCount++;
                    else {
                        // Look for the closest representative.
                        RepGenomeDb.Representation result = this.repDb.findClosest(seedFound);
                        String repId;
                        String repName;
                        String outlierFlag = "";
                        int sim = 0;
                        double dist = 1.0;
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
                            sim = result.getSimilarity();
                            dist = result.getDistance();
                            if (sim < repDb.getThreshold()) {
                                outCount++;
                                outlierFlag = "*";
                            }
                        }
                        if (! outlierFlag.isEmpty()) {
                            // Here we have an outlier.
                            if (this.updateFile != null) {
                                // Here we are updating the database and we need to add this outlier.
                                repName = genome.getName();
                                RepGenome rep = new RepGenome(seedFid, repName, seedFound);
                                this.repDb.addRep(rep);
                                repId = genome.getId();
                                sim = 9999;
                                dist = 0.0;
                                log.info("{} added to representative genome database.", genome);
                            }
                            if (this.saveDir != null) {
                                // Here we are saving the outliers.
                                this.target.add(genome);
                                log.info("{} saved to {}.", genome, this.saveDir);
                            }
                        }
                        // Write the output line.
                        writer.format("%s\t%s\t%s\t%s\t%d\t%8.4f\t%s%n", genome.getId(), genome.getName(), repId, repName,
                                sim, dist, outlierFlag);
                    }
                }
            }
            log.info("{} outliers found, {} skipped due to missing seed protein.",
                    outCount, skipCount);
            if (this.updateFile != null) {
                log.info("Writing new database to {}.", this.updateFile);
                repDb.save(this.updateFile);
            }
        } finally {
            if (this.outFile != null)
                this.outStream.close();
        }
    }

}
