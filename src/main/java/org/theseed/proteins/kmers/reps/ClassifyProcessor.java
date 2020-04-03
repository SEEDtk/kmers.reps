/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.kmers.KmerCollectionGroup;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;
import org.theseed.utils.ICommand;

/**
 * This processor builds an internal database from multiple FASTA files and then, for each input sequence
 * determines how far it is from a particular FASTA file in the database.
 *
 * The standard input is the FASTA file of sequences to process again the original files.  The comment
 * for each sequence is the list name of the FASTA file to compare it to.
 *
 * A control file contains the list name and file name of each original FASTA file.  This file should be
 * tab-delimited with no headers, the list name in the first column and the file name in the second.
 *
 * The output is tab-delimited, and consists of the input sequence ID, the closest-list name, and the distance.
 *
 * The positional parameter is the name of the control file.
 *
 * The command-line parameters are as follows.
 *
 * -K		protein kmer size (default 8)
 * -v		display progress messages on STDERR
 * -i		input FASTA file (default is STDIN)
 * -m		maximum acceptable distance
 *
 * @author Bruce Parrello
 *
 */
public class ClassifyProcessor implements ICommand {

    // FIELDS

    /** input FASTA stream */
    FastaInputStream inStream;
    /** input control file stream */
    TabbedLineReader controlStream;
    /** list of sequence lists */
    KmerCollectionGroup groups;
    /** output FASTA stream */
    FastaOutputStream outStream;

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

    /** maximum acceptable distance */
    @Option(name = "-m", aliases = { "--max", "--maxDist" }, usage = "maximum distance for a group")
    private double maxDist;

    /** control file */
    @Argument(index = 0, usage = "control file name", required = true)
    private File controlFile;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        this.inFile = null;
        this.maxDist = 0.995;
        this.setKmer(8);
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Set up the input FASTA.
                if (this.inFile != null)
                    this.inStream = new FastaInputStream(this.inFile);
                else
                    this.inStream = new FastaInputStream(System.in);
            }
            // Open the control file.
            this.controlStream = new TabbedLineReader(this.controlFile, 2);
            // Create the group controller.
            this.groups = new KmerCollectionGroup();
            // Denote we can proceed.
            retVal = true;
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        } catch (IOException e) {
            e.printStackTrace(System.err);
        }
        return retVal;
    }

    @Override
    public void run() {
        try {
            if (this.debug) System.err.println("Reading control file.");
            for (TabbedLineReader.Line line : this.controlStream) {
                String grp = line.get(0);
                File fastaFile = new File(line.get(1));
                if (this.debug) System.err.println("Reading group " + grp + " from " +
                        fastaFile + ".");
                try (FastaInputStream fastaStream = new FastaInputStream(fastaFile)) {
                    for (Sequence seq : fastaStream)
                        this.groups.addSequence(seq, grp);
                }
            }
            // Write the output header.
            System.out.println("seq_id\tgroup\tdistance");
            // Process the input file.  Note we only output distances less than the max.
            if (this.debug) System.err.println("Processing input sequences.");
            for (Sequence seq : this.inStream) {
                double distance = this.groups.getDistance(seq, seq.getComment());
                if (distance <= this.maxDist) {
                    System.out.format("%s\t%s\t%12.8f%n", seq.getLabel(), seq.getComment(),
                            distance);
                }
            }
        } catch (IOException e) {
            e.printStackTrace(System.err);
        } finally {
            this.controlStream.close();
            this.inStream.close();
        }
    }

}
