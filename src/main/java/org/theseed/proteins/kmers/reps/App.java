package org.theseed.proteins.kmers.reps;

import java.util.Arrays;

import org.theseed.genome.MD5Processor;
import org.theseed.utils.ICommand;

/**
 * This program processes protein kmers.  The commands are as follows.
 *
 * 	repdb		Create a representative-genome database from a FASTA file of protein sequences.
 *	group		Analyze proteins and group them together.
 *  classify	Compare proteins to multiple protein lists
 *  roles		Process a universal-role file against a representative genome database
 *  taxon		Create taxonomic training sets for each group
 *  genomes		Process genome evaluation results
 *  update		Process incremental genome evaluation results
 *  md5			Compute genome MD5s
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
        case "genomes" :
            processor = new GenomeProcessor();
            break;
        case "update" :
            processor = new UpdateProcessor();
            break;
        case "md5" :
            processor = new MD5Processor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
