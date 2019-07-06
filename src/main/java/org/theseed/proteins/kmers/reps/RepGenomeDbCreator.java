package org.theseed.proteins.kmers.reps;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

public class RepGenomeDbCreator {

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** similarity threshold */
    @Option(name="-m", aliases={"--sim", "--minScore"}, metaVar="100", usage="similarity threshold for representation")
    private int threshold;

    /** kmer size */
    @Option(name="-K", aliases={"--kmer"}, metaVar="8", usage="protein kmer size")
    private int kmerSize;
    
    
    public boolean parseCommand(String[] args) {
        // TODO Auto-generated method stub
        return false;
    }

    public void run() {
        // TODO Auto-generated method stub

    }

}
