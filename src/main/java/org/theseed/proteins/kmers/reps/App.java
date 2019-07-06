package org.theseed.proteins.kmers.reps;

/**
 * Create a representative-genome database from a FASTA file of protein sequences.
 *
 */
public class App
{
    public static void main( String[] args )
    {
        RepGenomeDbProcessor runObject = new RepGenomeDbProcessor();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
