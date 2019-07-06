package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * This is the primary class for creating a representative-genome database.  A FASTA file
 * is read from the standard input:  each record must have the key protein's feature ID as
 * the label, the genome name as the comment, and the key protein's sequence as the sequence.
 *
 * The positional parameter is the name of the output file.  This will be a FASTA file with
 * an empty-sequence header record, and one record per representative genome.  The command-line
 * options are as follows.
 *
 * -m	minimum similarity threshold for representation (default 100)
 * -K	protein kmer size (default 8)
 * -p	key protein name (default "Phenylalanyl-tRNA synthetase alpha chain")
 *
 * @author Bruce Parrello
 *
 */
public class RepGenomeDbCreator {

    // FIELDS
    /** number of genomes processed */
    private int genomesProcessed;
    /** start time of run */
    private long startTime;
    /** database being created */
    private RepGenomeDb repDB;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** similarity threshold */
    @Option(name="-m", aliases={"--sim", "--minScore"}, metaVar="100", usage="similarity threshold for representation")
    private int threshold;

    /** kmer size */
    @Option(name="-K", aliases={"--kmer"}, metaVar="8", usage="protein kmer size")
    private void setKmer(int newSize) {
        ProteinKmers.setKmerSize(newSize);
    }

    /** key protein name */
    @Option(name="-p", aliases={"--prot", "--role"}, metaVar="\"role name\"", usage="name of the role for the key protein")
    private String protName;

    /** output file */
    @Argument(index=0, metaVar="outFile", usage="output result file",
            required=true, multiValued=false)
    private File outFile;

    /** Parse the command line parameters and options. */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.threshold = 100;
        this.protName = RepGenomeDb.DEFAULT_PROTEIN;
        this.outFile = null;
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

    /** Create the representative-genome database. */
    public void run() {
        try {
            System.out.format("Creating rep db with K=%d and threshold %d using %s.%n",
                    ProteinKmers.kmerSize(), this.threshold, this.protName);
            repDB = new RepGenomeDb(this.threshold, this.protName);
            FastaInputStream inStream = new FastaInputStream(System.in);
            // Loop through the input stream, processing representative genomes.  We can
            // do this in one statement, but we slow down to display progress.
            System.out.println("Processing input.");
            genomesProcessed = 0;
            long lastTime = System.currentTimeMillis();
            startTime = lastTime;
            for (Sequence inSequence : inStream) {
                RepGenome newGenome = new RepGenome(inSequence);
                repDB.checkGenome(newGenome);
                genomesProcessed++;
                if ((System.currentTimeMillis() - lastTime) > 10000) {
                    // A minute since our last progress. Show more progress.
                    showProgress();
                    lastTime = System.currentTimeMillis();
                }
            }
            showProgress();
            inStream.close();
            // Save the database.
            System.out.println("Saving database to " + outFile.getPath() + ".");
            repDB.save(outFile);
        } catch (IOException e) {
            throw new RuntimeException("Error reading FASTA file.", e);
        }
    }

    private void showProgress() {
        long rate = 0;
        if (genomesProcessed > 0) {
            rate = genomesProcessed * 1000 / (System.currentTimeMillis() - startTime);
        }
        System.out.format("%d total genomes input, %d kept (%d genomes/second).%n",
                genomesProcessed, repDB.size(), rate);
    }
}

