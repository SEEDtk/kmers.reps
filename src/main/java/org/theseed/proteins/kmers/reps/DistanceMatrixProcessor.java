/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.utils.ICommand;

/**
 * This class reads a representative-genome database and computes the distance between each pair of genomes
 * in the database.  The positional parameter is the name of the representative-genome database file.  The
 * output will be a three-column table on the standard output.
 *
 * The command-line parameters are as follows.
 *
 * -h 	display command-line usage
 * -v	display progress messages on STDERR
 *
 * @author Bruce Parrello
 *
 */
public class DistanceMatrixProcessor implements ICommand {

    // FIELDS
    /** target representative-genome database */
    private RepGenomeDb repDB;

    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** produce progress messages on STDERR */
    @Option(name="--verbose", aliases={"-v", "--debug"}, usage="show progress on STDERR")
    private boolean debug;

    /** representative-genome file */
    @Argument(index=0, metaVar="repDbFile", usage="representative genome database file",
            required=true, multiValued=false)
    private File dbFile;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        this.help = false;
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else if (! this.dbFile.canRead()) {
                throw new FileNotFoundException("RepDB file " + this.dbFile + " not found or unreadable.");
            } else {
                retVal = true;
            }
        } catch (CmdLineException | FileNotFoundException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        }
        return retVal;
    }

    @Override
    public void run() {
        try {
            // Load the rep-genome database.
            if (debug) System.err.println("Loading representative genome data from " + this.dbFile.getPath());
            this.repDB = RepGenomeDb.load(this.dbFile);
            // Write the output header.
            System.out.println("genome_1\tgenome_2\tdistance");
            // Compute the number of pairs to process.
            int totalPairs = this.repDB.size() * (this.repDB.size() + 1) / 2;
            if (debug) System.err.format("%d total genomes, requiring %d pairs to be processed.%n", this.repDB.size(), totalPairs);
            // Loop through the representatives.
            long startTime = System.currentTimeMillis();
            int processed = 0;
            for (RepGenome rep : this.repDB) {
                for (RepGenome rep2 : this.repDB) {
                    // Compare the two genomes to insure each pair is only processed one.
                    if (rep.compareTo(rep2) < 0) {
                        double distance = rep.distance(rep2);
                        System.out.format("%s\t%s\t%8.6f%n", rep.getGenomeId(), rep2.getGenomeId(), distance);
                        processed++;
                        if (debug && processed % 5000 == 0) {
                            double rate = (double) processed / ((System.currentTimeMillis() - startTime) / 1000.0);
                            System.err.format("%d of %d pairs processed, %5.1f genomes/second.%n", processed, totalPairs,
                                    rate);
                        }
                    }
                }
            }
            long duration = System.currentTimeMillis() - startTime;
            System.err.format("%d pairs processed in %d seconds.", processed, duration);
        } catch (IOException e) {
            throw new RuntimeException("Error reading FASTA file.", e);
        }
    }

}
