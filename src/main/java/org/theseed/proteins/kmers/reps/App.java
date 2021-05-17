package org.theseed.proteins.kmers.reps;

import java.util.Arrays;

import org.theseed.proteins.UniRoleProcessor;
import org.theseed.genome.MD5Processor;
import org.theseed.genome.coupling.CouplesProcessor;
import org.theseed.genome.coupling.PrepareProcessor;
import org.theseed.utils.BaseProcessor;

/**
 * This program processes protein kmers.  The commands are as follows.
 *
 *  gtoReps		find the representative genome of each GTO in a directory
 * 	repdb		Create a representative-genome database from a FASTA file of protein sequences.
 *	group		Analyze proteins and group them together.
 *  classify	Compare proteins to multiple protein lists
 *  genomes		Process genome evaluation results
 *  update		Process incremental genome evaluation results
 *  md5			Compute genome MD5s
 *  coupling	Compute functional coupling for a set of genomes
 *  prepare		Prepare a GTO directory for use in the coupling website
 *  distances	Create a distance matrix for representative genomes
 *  seqTable	Create a table of identifying sequences for each representative genome
 *  seqComp		Create a table comparing PheS distance to SSU-rRNA distance
 *  seqTest		Compare the closest PheS genome to the closest SSU-rRNA genome
 *  target		Find a kmer target in a set of genomes
 *  univ		create a report on the singly-occurring roles in a group of genomes
 *  gtoClass	Find representatives for GTOs in multiple RepGen databases
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        switch (command) {
        case "gtoReps" :
            processor = new GtoRepGenomeProcessor();
            break;
        case "distances" :
            processor = new DistanceMatrixProcessor();
            break;
        case "repdb" :
            processor = new RepGenomeDbProcessor();
            break;
        case "group" :
            processor = new RepMatrixProcessor();
            break;
        case "classify" :
            processor = new ClassifyProcessor();
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
        case "coupling" :
            processor = new CouplesProcessor();
            break;
        case "prepare" :
            processor = new PrepareProcessor();
            break;
        case "seqTable" :
            processor = new SeqTableProcessor();
            break;
        case "seqComp" :
            processor = new SeqCompProcessor();
            break;
        case "seqTest" :
            processor = new SeqTestProcessor();
            break;
        case "target" :
            processor = new TargetProcessor();
            break;
        case "univ" :
            processor = new UniRoleProcessor();
            break;
        case "gtoClass" :
            processor = new GtoClassProcessor();
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
