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
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.genome.iterator.GenomeTargetType;
import org.theseed.genome.iterator.IGenomeTarget;
import org.theseed.io.TabbedLineReader;
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
 * --filter		if specified, a tab-delimited file of genome IDs; only the genome whose IDs are listed in the
 * 				first column will be processed, and they will be processed in the order presented
 * --blackList	if specified, a tab-delimited file of genome IDs; the genomes whose IDs are specified will
 * 				be processed last
 * --rna		if specified, genomes without a valid SSU rRNA will be skipped; genomes with a short SSU rRNA
 * 				will be processed before the blacklist
 * --update		if specified, a file name; the repGenomeDB will be updated and stored in the named file
 * --save		if specified, a genome target to which outliers should be saved
 * --saveType	type of genome target (MASTER or DIR)
 * --create		create a new representative-genome database with default parameters and the specified similarity threshold
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
    /** input genome source */
    private GenomeSource source;
    /** genome ID filter list */
    private List<String> idList;
    /** output genome target */
    private IGenomeTarget target;
    /** output stream */
    private OutputStream outStream;
    /** list of genome IDs to defer */
    private Set<String> blacklist;
    /** deferred-genome list */
    private List<RepGenome> deferrals;
    /** outlier count */
    private int outCount;
    /** minimum acceptable rRNA length */
    private static final int MIN_RNA_LEN = 1400;

    // COMMAND-LINE OPTIONS

    /** similarity threshold; indicates a new database will be created */
    @Option(name = "--create", metaVar = "200", usage = "if specified, indicates that the rep-genome database should be created with the indicated threshold")
    private int createFlag;

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

    /** blacklist file */
    @Option(name = "--blackList", metaVar = "blacklistFile.tbl", usage = "if specified, a file of genome IDs to defer")
    private File blackListFile;

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
        this.createFlag = 0;
        this.blackListFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (! this.gtoDir.exists())
            throw new FileNotFoundException("Input source " + this.gtoDir + " not found.");
        // Load the genome source.
        this.source = this.type.create(this.gtoDir);
        log.info("{} genomes in source {}.", this.source.size(), this.gtoDir);
        // Load the rep-genome database.
        if (this.createFlag == 0) {
            log.info("Loading representative-genome database from {}.", this.repDbFile);
            this.repDb = RepGenomeDb.load(this.repDbFile);
            // Report on the update status.
            if (this.updateFile != null)
                log.info("Outliers will be added to the database and the updated database written to {}.", this.updateFile);
        } else if (this.updateFile != null)
            throw new ParseFailureException("Update file is mutually exclusive with create mode.");
        else {
            // Set up to write the resulting database to the specified rep-genome file.
            log.info("Creating new representative-genome database in {}.", this.repDbFile);
            this.repDb = new RepGenomeDb(this.createFlag);
            this.updateFile = this.repDbFile;
        }
        // If there is a filter file, create the filter set.
        if (this.filterFile == null) {
            log.info("All input genomes will be processed.");
            this.idList = new ArrayList<String>(this.source.getIDs());
        } else {
            this.idList = TabbedLineReader.readColumn(this.filterFile, "1");
            log.info("{} genomes specified in filter file {}.", this.idList.size(), this.filterFile);
        }
        // If there is a blacklist file, set up the black list.
        if (this.blackListFile == null) {
            this.blacklist = Collections.emptySet();
        } else {
            this.blacklist = TabbedLineReader.readSet(this.blackListFile, "1");
            log.info("{} genomes in blacklist read from {}.", this.blacklist.size(), this.blackListFile);
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
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        try (PrintWriter writer = new PrintWriter(this.outStream)) {
            // Start the output report.
            writer.println("genome_id\tgenome_name\trep_id\trep_name\tsimilarity\tdistance\toutlier");
            // Set up the deferral list.
            this.deferrals = new ArrayList<RepGenome>(2000);
            // Now read through the genome directory.
            log.info("Scanning genome directory {}.", this.gtoDir);
            int skipCount = 0;
            this.outCount = 0;
            int gCount = 0;
            for (String genomeId : this.idList) {
                gCount++;
                Genome genome = this.source.getGenome(genomeId);
                // It is completely acceptable for the genome to not be in the source.
                // We only process if it's found.
                if (genome != null) {
                    log.info("Processing genome {} of {}: {}.", gCount, this.idList.size(), genome);
                    RepGenome rep = this.repDb.getSeedProtein(genome);
                    // Only proceed if we have a seed protein.  Here we also check for the SSU rRNA.
                    boolean skip = false;
                    if (rep == null) {
                        log.warn("{} does not have a seed protein.", genome);
                        skip = true;
                        skipCount++;
                    } else if (this.blacklist.contains(genome.getId())) {
                        log.info("{} deferred due to blacklist.", genome);
                        this.defer(rep);
                        skip = true;
                    } else if (this.rnaFlag) {
                        String ssuRRNA = genome.getSsuRRna();
                        if (ssuRRNA.isEmpty()) {
                            log.warn("{} does not have an SSU rRNA.", genome);
                            skip = true;
                            skipCount++;
                        } else if (ssuRRNA.length() < MIN_RNA_LEN) {
                            log.warn("{} deferred due to short RNA.", genome);
                            this.defer(rep);
                            skip = true;
                        }
                    }
                    if (! skip) {
                        // Look for the closest representative.
                        processGenome(writer, rep);
                    }
                }
            }
            // Now process the deferrals.
            log.info("Processing {} deferred genomes.", this.deferrals.size());
            for (RepGenome deferral : this.deferrals) {
                processGenome(writer, deferral);
            }
            log.info("{} outliers found, {} skipped due to missing seed protein.",
                    this.outCount, skipCount);
            if (this.updateFile != null) {
                log.info("Writing new database to {}.", this.updateFile);
                repDb.save(this.updateFile);
            }
            // Insure the output directory (if any) is cleanly finished.
            if (this.target != null)
                this.target.finish();
        } finally {
            if (this.outFile != null)
                this.outStream.close();
        }
    }

    /**
     * Process a genome to find its representative.
     *
     * @param writer		print writer for the output line
     * @param rep			representation object for genome to process
     *
     * @throws IOException
     */
    protected void processGenome(PrintWriter writer, RepGenome rep) throws IOException {
        RepGenomeDb.Representation result = this.repDb.findClosest(rep);
        String repId;
        String repName;
        String outlierFlag = "";
        int sim = 0;
        double dist = 1.0;
        if (result.getSimilarity() == 0) {
            // Here there is nothing close.
            repId = "";
            repName = "";
            this.outCount++;
            outlierFlag = "*";
        } else {
            // Compute the ID and name of the closest rep.
            repId = result.getGenomeId();
            repName = result.getRepresentative().getName();
            sim = result.getSimilarity();
            dist = result.getDistance();
            if (sim < repDb.getThreshold()) {
                this.outCount++;
                outlierFlag = "*";
            }
        }
        if (! outlierFlag.isEmpty()) {
            // Here we have an outlier.
            if (this.updateFile != null) {
                // Here we are updating the database and we need to add this outlier.
                this.repDb.addRep(rep);
                repId = rep.getGenomeId();
                sim = 9999;
                dist = 0.0;
                log.info("{} ({}) added to representative genome database.", rep.getGenomeId(), rep.getName());
            }
            if (this.saveDir != null) {
                // Here we are saving the outliers.
                Genome genome = this.source.getGenome(rep.getGenomeId());
                this.target.add(genome);
                log.info("{} saved to {}.", genome, this.saveDir);
            }
        }
        // Write the output line.
        writer.format("%s\t%s\t%s\t%s\t%d\t%8.4f\t%s%n", rep.getGenomeId(), rep.getName(), repId, repName,
                sim, dist, outlierFlag);
    }

    /**
     * Defer a genome for later processing.
     *
     * @param genome		genome of interest
     * @param seedFid		seed protein feature ID
     * @param seedFound		seed protein sequence
     */
    private void defer(RepGenome rep) {
        this.deferrals.add(rep);
    }

}
