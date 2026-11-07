package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.seeds.FinderKmerRepGenerator;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This processor generates representative sets for all the genome groups listed in a directory.
 * The input is a directory containing genome group files, and the output is a set of representative sets for each group.
 * Each genome group file should have the extension ".tbl", and it must be tab-delimited with headers. The genome ID is
 * generally taken from the first column.
 * 
 * The output for each genome group will be a representative set file, typically with the same base name as the input file,
 * but the modifier "repgen" inserted before the file extension. All the output files will be put in the output directory.
 * 
 * The positional parameter is the name of the input directory containing the genome group files.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -D   output directory for the representative set files (default "RepGen" in the current directory)
 * -c   index (1-based) or name of input column containing the genome IDs
 * -N   target size for the output repgen set (default 200)
 * -b   batch size for processing genome IDs (default 400)
 * 
 * --roles  role definition file for the roles to use (default roles.for.finder in the current directory)
 * --fuzz   a distance from the computed closeness limit, relative to the standard deviation, for alternate limits to try (default is 0.3)
 *          A fuzz factor of 0 means the fuzz factor will not be applied.
* 
 * MassRepFinderProcessor
 */
public class MassRepFinderProcessor extends BaseMultiReportProcessor {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(MassRepFinderProcessor.class);
    /** generator for creating the repgen sets */
    private FinderKmerRepGenerator generator;
    /** input files to read */
    private Collection<File> inputFiles;
    /** file extension to use for the genome files */
    private static final String[] GENOME_FILE_EXTENSIONS = new String[] { "tbl" };

    /** index (1-based) or name of input column containing the genome IDs */
    @Option(name = "--col", aliases = { "-c" }, metaVar = "genome_id", usage = "index (1-based) or name of input column containing the genome IDs")
    private String genomeIdCol;

    /** distance from the computed closeness limit, relative to the standard deviation, for alternate limits to try */
    @Option(name = "--fuzz", metaVar = "0.2", usage = "distance from the computed closeness limit, relative to the standard deviation, for alternate limits to try")
    private double fuzz;

    /** role definition file for the roles to use */
    @Option(name = "--roles", metaVar = "roles.for.finder", usage = "role definition file for the roles to use")
    private File rolesFile;

    /** batch size for processing genome IDs */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "300", usage = "batch size for processing genome IDs")
    private int batchSize;

    /** target number of genomes for the output representative-genome set */
    @Option(name = "--target", aliases = { "-N" }, metaVar = "1000", usage = "target number of genomes for the output representative-genome set")
    private int targetSize;

    /** input directory for genome files */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory for genome files")
    private File inDir;


    @Override
    protected File setDefaultOutputDir(File curDir) {
        return new File(curDir, "RepGens");
    }

    @Override
    protected void setMultiReportDefaults() {
        this.fuzz = 0.3;
        this.genomeIdCol = "1";
        this.rolesFile = new File("roles.for.finder");
        this.targetSize = 200;
        this.batchSize = 400;
    }
    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        // Validate the target size.
        if (this.targetSize <= 0)
            throw new ParseFailureException("Target size must be positive.");
        // Validate the fuzz factor.
        if (this.fuzz < 0.0)
            throw new ParseFailureException("Fuzz factor cannot be negative.");
        // Insure the role map is readable.
        if (! this.rolesFile.canRead())
            throw new ParseFailureException("Role definition file " + this.rolesFile + " is not found or unreadable.");
        // Insure the input directory is readable.
        if (! this.inDir.isDirectory() || ! this.inDir.canRead())
            throw new IOException("Input directory " + this.inDir + " is not found or unreadable.");
        // Scan the input directory for genome files.
        this.inputFiles = FileUtils.listFiles(this.inDir, GENOME_FILE_EXTENSIONS, false);
        if (this.inputFiles.isEmpty())
            throw new IOException("No genome files found in input directory " + this.inDir + ".");
        log.info("{} genome files found in {}.", this.inputFiles.size(), this.inDir);
    }
    @Override
    protected void runMultiReports() throws Exception {
        // Create the repgen generator.
        this.generator = new FinderKmerRepGenerator(this.rolesFile, this.fuzz, this.batchSize, this.targetSize);
        // Loop through the input files, processing each group.
        int fileCount = 0;
        for (File inFile : this.inputFiles) {
            fileCount++;
            log.info("Processing genome file {} of {}: {}.", fileCount, this.inputFiles.size(), inFile);
            String inFileName = inFile.getName();
            String outFileName = StringUtils.substringBeforeLast(inFileName, ".") + ".repgen.tbl";
            File outFile = this.getOutFile(outFileName);
            try (TabbedLineReader inStream = new TabbedLineReader(inFile);
                 PrintWriter outStream = new PrintWriter(outFile)) {
                this.generator.processGenomes(inStream, this.genomeIdCol, outStream);
            }
        }
        log.info("Processing complete. {} genome files processed.", fileCount);
    }

}
