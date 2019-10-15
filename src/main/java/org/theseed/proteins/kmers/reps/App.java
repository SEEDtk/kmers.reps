package org.theseed.proteins.kmers.reps;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * This program processes protein kmers.  The commands are as follows.
 *
 * 	repdb		Create a representative-genome database from a FASTA file of protein sequences.
 *	group		Analyze proteins and group them together.
 *  classify	Compare proteins to multiple protein lists
 *  roles		Process a universal-role file against a representative genome database
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        switch (command) {
        case "repdb" :
            processor = new RepGenomeDbProcessor();
            break;
        case "group" :
            processor = new RepMatrixProcessor();
            break;
        case "classify" :
            processor = new ClassifyProcessor();
            break;
        case "roles" :
            processor = new RolesProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ": must be \"repdb\" or \"matrix\".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
