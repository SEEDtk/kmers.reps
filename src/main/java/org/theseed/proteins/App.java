package org.theseed.proteins;

import org.theseed.proteins.kmers.reps.RepGenomeDbCreator;

/**
 * Create a representative-genome database from a FASTA file of protein sequences.
 *
 */
public class App
{
    public static void main( String[] args )
    {
        RepGenomeDbCreator runObject = new RepGenomeDbCreator();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
