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
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.reports.RepTaxonReport;

/**
 * This command reads a genome source and outputs the ID, name, and high-level taxonomy for each genome.  The taxonomy can be
 * tweaked to start with species, genus, or family.
 *
 * The positional parameter is the file or directory name for the genome source.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more detailed log messages
 * -o	output file for report (if not STDOUT)
 * -t	type of genome source (default DIR)
 *
 * --rank	minimum rank to display (default SPECIES)
 *
 * @author Bruce Parrello
 *
 */
public class RepTaxProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(RepTaxProcessor.class);
    /** input genome source */
    private GenomeSource genomes;

    // COMMAND-LINE OPTIONS

    /** genome source type */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** minimum display rank */
    @Option(name = "--rank", usage = "minimum rank to display")
    private RepTaxonReport.MinRank minRank;

    /** genome source */
    @Argument(index = 0, metaVar = "inDir", usage = "genome source directory or file", required = true)
    private File inDir;


    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.minRank = RepTaxonReport.MinRank.SPECIES;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " is not found.");
        // Connect to the genome source.
        this.genomes = this.sourceType.create(this.inDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.inDir);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Start a timer for tracing.
        long start = System.currentTimeMillis();
        int gCount = 0;
        // Create the report helper.
        var reporter = new RepTaxonReport(writer, this.minRank);
        for (Genome genome : this.genomes) {
            reporter.write(genome);
            gCount++;
            if (log.isInfoEnabled() && System.currentTimeMillis() - start > 5000) {
                log.info("{} of {} genomes output.", gCount, this.genomes.size());
                start = System.currentTimeMillis();
            }
        }
    }

}
