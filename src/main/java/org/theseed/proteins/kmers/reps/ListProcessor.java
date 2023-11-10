/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.utils.BaseReportProcessor;

/**
 * This is a simple command that lists the ID and name of each genome in a representative-genome database.
 *
 * The positional parameter is the name of the representative-genome database.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	name of the output file (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class ListProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ListProcessor.class);
    /** repgen database */
    private RepGenomeDb repDb;

    // COMMAND-LINE OPTIONS

    /** repgen DB file name */
    @Argument(index = 0, metaVar = "repDb.ser", usage = "name of the representative-genome database file", required = true)
    private File rebDbFile;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify the repgen database.
        if (! this.rebDbFile.canRead())
            throw new FileNotFoundException("Repgen DB file " + this.rebDbFile + " is not found or unreadable.");
        log.info("Loading repgen database from {}.", this.rebDbFile);
        this.repDb = RepGenomeDb.load(this.rebDbFile);
        log.info("{} representatives loaded.", this.repDb.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Sort the genomes.
        RepGenome[] reps = this.repDb.all();
        // Write the header line.
        writer.println("genome_id\tgenome_name");
        // Loop through the genomes.
        for (RepGenome rep : reps)
            writer.println(rep.getGenomeId() + "\t" + rep.getName());
    }

}
