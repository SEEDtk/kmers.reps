/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This script takes two files containing genome ID lists and sorts the second file so the genome IDs in the first file
 * occur first.  The output file will contain only the genome IDs.
 *
 * The positional parameters are the names of the two files-- the priority file and the main file.  The command-line
 * options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file name (if not STDOUT)
 * -c	column ID for genome IDs in the second file (default "1")
 *
 * --col1	column ID for genome IDs in the first file (default "1")
 *
 * @author Bruce Parrello
 *
 */
public class PrioritizeProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PrioritizeProcessor.class);
    /** set of IDs in the first file */
    private Set<String> priorityIDs;
    /** list of IDs in the second file that are not in the first file */
    private List<String> secondIDs;
    /** list of IDs in the first file that are in the second file */
    private List<String> firstIDs;

    // COMMAND-LINE OPTIONS

    /** column ID for the first (priority) file */
    @Option(name = "--col1", metaVar = "rep_id", usage = "index (1-based) or name of ID column in first file")
    private String col1Spec;

    /** column ID for the second (list) file */
    @Option(name = "--col2", aliases = { "-c" }, metaVar = "genome_id", usage = "index (1-based) or name of ID column in second file")
    private String col2Spec;

    /** priority input file */
    @Argument(index = 0, metaVar = "priorityFile.tbl", usage = "file containing priority IDs")
    private File priorityFile;

    /** list input file */
    @Argument(index = 1, metaVar = "listFile.tbl", usage = " file containing full list of IDs")
    private File listFile;

    @Override
    protected void setReporterDefaults() {
        this.col1Spec = "1";
        this.col2Spec = "1";
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify both input files.
        if (! this.priorityFile.canRead())
            throw new FileNotFoundException("Priority input file " + this.priorityFile + " is not found or unreadable.");
        if (! this.listFile.canRead())
            throw new FileNotFoundException("List input file " + this.listFile + " is not found or unreadable.");
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Read in the priority IDs.
        this.priorityIDs = TabbedLineReader.readSet(this.priorityFile, this.col1Spec);
        log.info("{} IDs read from priority file.", this.priorityIDs.size());
        // We need to loop through the list file, putting the IDs found into two lists.
        this.firstIDs = new ArrayList<String>(this.priorityIDs.size());
        this.secondIDs = new ArrayList<String>((int) (this.listFile.length() / 50) + this.priorityIDs.size());
        try (TabbedLineReader inStream = new TabbedLineReader(this.listFile)) {
            // Get the index and name of the ID column.
            int colIdx = inStream.findField(this.col2Spec);
            String colName = inStream.getLabels()[colIdx];
            // The output file has just the ID column as its header.
            writer.println(colName);
            // Now read all the IDs, in order.
            for (TabbedLineReader.Line line : inStream) {
                String id = line.get(colIdx);
                if (this.priorityIDs.contains(id))
                    this.firstIDs.add(id);
                else
                    this.secondIDs.add(id);
            }
        }
        log.info("{} IDs read from list file.", this.firstIDs.size() + this.secondIDs.size());
        // Now unspool the results.
        for (String id : this.firstIDs)
            writer.println(id);
        for (String id : this.secondIDs)
            writer.println(id);
    }


}
