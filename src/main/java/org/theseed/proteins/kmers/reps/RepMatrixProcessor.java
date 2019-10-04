/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.utils.ICommand;

/**
 * This command reads a FASTA file and creates an internal distance matrix for the sequences.  The command will compute
 * the centermost sequence and output the minimum, mean, and maximum distance from that sequence.
 *
 * The parameters are as follows.
 *
 * -K		protein kmer size (default 8)
 * -v		display progress messages on STDERR
 * -i		input FASTA file (default is STDIN)
 *
 * @author Bruce Parrello
 *
 */
public class RepMatrixProcessor implements ICommand {

    // FIELDS

    /** input FASTA stream */
    FastaInputStream inStream;

    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name = "-h", aliases = { "--help" }, help = true)
    protected boolean help;

    /** TRUE if we want progress messages */
    @Option(name = "-v", aliases = { "--verbose", "--debug" }, usage = "display progress on STDERR")
    protected boolean debug;

    /** kmer size */
    @Option(name="-K", aliases={"--kmer"}, metaVar="8", usage="protein kmer size (default 8)")
    private void setKmer(int newSize) {
        ProteinKmers.setKmerSize(newSize);
    }

    /** input file (if not using STDIN) */
    @Option(name = "--input", aliases = { "-i" }, usage = "input file (if not STDIN)")
    private File inFile;


    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        this.inFile = null;
        this.setKmer(8);
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {

                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        }
        return retVal;
    }

    @Override
    public void run() {
        // TODO read input, create list of sequenceInfos
    }

}
