package org.theseed.proteins.kmers.reps;

import java.util.Arrays;

import org.theseed.proteins.UniRoleProcessor;
import org.theseed.genome.MD5Processor;
import org.theseed.genome.coupling.CouplesProcessor;
import org.theseed.genome.coupling.PrepareProcessor;
import org.theseed.sequence.RnaVerifyProcessor;
import org.theseed.utils.BaseProcessor;

/**
 * This program processes protein kmers.  The commands are as follows.
 *
 *  build		Create a representative-genome database for a pre-selected set of representatives
 *  gtoReps		find the representative genome of each GTO in a directory
 *  ssuReps		find the representative genome of each GTO in a directory using an SSU table
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
 *  ssuCheck	Find bad SSU rRNA sequences by blasting against SILVA
 *  target		Find a kmer target in a set of genomes
 *  univ		create a report on the singly-occurring roles in a group of genomes
 *  gtoClass	Find representatives for GTOs in multiple RepGen databases
 *  missing		Produce the missing-roles report for a directory of GTOs
 *  prio		prioritize a list of genomes using a second list
 *  list		list the genomes in a representative genome database
 *  rnaVerify	build a blacklist of genomes with bad or suspicious SSU rRNA sequences
 *  seedTable	build a 4-column table for a specified seed protein
 *  buildRefDb	build a reference-genome database for evaluation
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
        case "ssuReps" :
            processor = new SsuRepGenomeProcessor();
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
        case "missing" :
            processor = new MissingRoleProcessor();
            break;
        case "gtoClass" :
            processor = new GtoClassProcessor();
            break;
        case "ssuCheck" :
            processor = new BadRnaProcessor();
            break;
        case "prio" :
            processor = new PrioritizeProcessor();
            break;
        case "list" :
            processor = new ListProcessor();
            break;
        case "rnaVerify" :
            processor = new RnaVerifyProcessor();
            break;
        case "seedTable" :
            processor = new SeedTableProcessor();
            break;
        case "build" :
            processor = new BuildRepDbProcessor();
            break;
        case "buildRefDb" :
            processor = new BuildRefDbProcessor();
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
