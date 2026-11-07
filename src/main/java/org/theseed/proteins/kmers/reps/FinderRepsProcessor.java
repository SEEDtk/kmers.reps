package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.kohsuke.args4j.Option;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.seeds.FinderKmerRepGenerator;
import org.theseed.utils.BasePipeProcessor;

/**
 * This subcommand generates a representative set of genomes using finder-kmer technology. This is a two-pass process. First,
 * we run through the genomes to get key metrics about the set. We estimate the closeness limit to use based on those
 * metrics. Then we make another pass to generate the representative set itself. A temporary file containing the protein
 * sequences discovered in the first pass will be used to speed up the second pass.
 * 
 * The genome list should be provided on the standard input as a tab-delimited file with headers. The default column for the
 * genome ID will be the first one. The representative-genome list will be written to the standard output.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -i   input file containing the genome list (if not STDIN)
 * -o   output file for the representative-genome list (if not STDOUT)
 * -c   index (1-based) or name of input column containing the genome IDs
 * -N   target size for the output repgen set (default 200)
 * -b   batch size for processing genome IDs (default 400)
 * 
 * --roles  role definition file for the roles to use (default roles.for.finder in the current directory)
 * --fuzz   a distance from the computed closeness limit, relative to the standard deviation, for alternate limits to try (default is 0.3)
 *          A fuzz factor of 0 means the fuzz factor will not be applied.
 */
public class FinderRepsProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(FinderRepsProcessor.class);
    /** generator for creating the repgen sets */
    private FinderKmerRepGenerator generator;
    
    // COMMAND-LINE OPTIONS

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
    
    @Override
    protected void setPipeDefaults() {
        this.fuzz = 0.3;
        this.genomeIdCol = "1";
        this.rolesFile = new File("roles.for.finder");
        this.targetSize = 200;
        this.batchSize = 400;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        int idx = inputStream.findField(this.genomeIdCol);
        log.info("Genome IDs will be taken from column {}.", idx + 1);
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Validate the target size.
        if (this.targetSize <= 0)
            throw new ParseFailureException("Target size must be positive.");
        // Validate the fuzz factor.
        if (this.fuzz < 0.0)
            throw new ParseFailureException("Fuzz factor cannot be negative.");
        // Insure the role map is readable.
        if (! this.rolesFile.canRead())
            throw new ParseFailureException("Role definition file " + this.rolesFile + " is not found or unreadable.");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        this.generator = new FinderKmerRepGenerator(this.rolesFile, this.fuzz, this.batchSize, this.targetSize);
        try {
            this.generator.processGenomes(inputStream, this.genomeIdCol, writer);
        } finally {
            log.info("Closing the generator.");
            this.generator.close();
        }
    }

}
