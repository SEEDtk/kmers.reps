package org.theseed.proteins.kmers.reps;

import java.util.Arrays;

import org.theseed.basic.BaseProcessor;
import org.theseed.genome.MD5Processor;
import org.theseed.genome.coupling.CouplesProcessor;
import org.theseed.genome.coupling.PrepareProcessor;
import org.theseed.proteins.UniRoleProcessor;
import org.theseed.sequence.RnaVerifyProcessor;

/**
 * This program processes protein kmers.  The commands are as follows.
 *
 *  build			Create a representative-genome database for a pre-selected set of representatives
 *  gtoReps			find the representative genome of each GTO in a directory
 *  ssuReps			find the representative genome of each GTO in a directory using an SSU table
 * 	repdb			Create a representative-genome database from a FASTA file of protein sequences.
 *	group			Analyze proteins and group them together.
 *  classify		Compare proteins to multiple protein lists
 *  genomes			Process genome evaluation results
 *  md5				Compute genome MD5s
 *  coupling		Compute functional coupling for a set of genomes
 *  prepare			Prepare a GTO directory for use in the coupling website
 *  distances		Create a distance matrix for representative genomes
 *  seqTable		Create a table of identifying sequences for each representative genome
 *  seqComp			Create a table comparing PheS distance to SSU-rRNA distance
 *  seqTest			Compare the closest PheS genome to the closest SSU-rRNA genome
 *  ssuCheck		Find bad SSU rRNA sequences by blasting against SILVA
 *  target			Find a kmer target in a set of genomes
 *  univ			create a report on the singly-occurring roles in a group of genomes
 *  gtoClass		Find representatives for GTOs in multiple RepGen databases
 *  missing			Produce the missing-roles report for a directory of GTOs
 *  prio			prioritize a list of genomes using a second list
 *  list			list the genomes in a representative genome database
 *  rnaVerify		build a blacklist of genomes with bad or suspicious SSU rRNA sequences
 *  seedTable		build a 4-column table for a specified seed protein
 *  buildRefDb		build a reference-genome database for evaluation
 *  updateMaster	update a master PATRIC database:  remove obsolete genomes and add the new ones
 *  taxReport		determine how much information a genus or species identification tells us about repgen membership
 *  repTax			write the taxonomy report for a genome source
 *  neighbors		list the N closest neighboring representative genomes to all genomes in a source
 *  distDist		find a set of genomes that are various distributed distances from a key genome
 *  buildFromGtos	build a representative genome database for all genomes in a source
 *  hiCheck         check the master representation list to see which genomes violate repgen hierarchy
 *  repSub          validate the repgen sets to insure they are hierarchical
 *  buildSets       create neighbor sets from a genome representation list file
 */
public class App
{
    /** static array containing command names and comments */
    protected static final String[] COMMANDS = new String[] {
             "build", "Create a representative-genome database for a pre-selected set of representatives",
             "gtoReps", "find the representative genome of each GTO in a directory",
             "ssuReps", "find the representative genome of each GTO in a directory using an SSU table",
             "repdb", "Create a representative-genome database from a FASTA file of protein sequences.",
             "group", "Analyze proteins and group them together.",
             "classify", "Compare proteins to multiple protein lists",
             "genomes", "Process genome evaluation results",
             "md5", "Compute genome MD5s",
             "coupling", "Compute functional coupling for a set of genomes",
             "prepare", "Prepare a GTO directory for use in the coupling website",
             "distances", "Create a distance matrix for representative genomes",
             "seqTable", "Create a table of identifying sequences for each representative genome",
             "seqComp", "Create a table comparing PheS distance to SSU-rRNA distance",
             "seqTest", "Compare the closest PheS genome to the closest SSU-rRNA genome",
             "ssuCheck", "Find bad SSU rRNA sequences by blasting against SILVA",
             "target", "Find a kmer target in a set of genomes",
             "univ", "create a report on the singly-occurring roles in a group of genomes",
             "gtoClass", "Find representatives for GTOs in multiple RepGen databases",
             "missing", "Produce the missing-roles report for a directory of GTOs",
             "prio", "prioritize a list of genomes using a second list",
             "list", "list the genomes in a representative genome database",
             "rnaVerify", "build a blacklist of genomes with bad or suspicious SSU rRNA sequences",
             "seedTable", "build a 4-column table for a specified seed protein",
             "buildRefDb", "build a reference-genome database for evaluation",
             "updateMaster", "update a master PATRIC database:  remove obsolete genomes and add the new ones",
             "taxReport", "determine how much information a genus or species identification tells us about repgen membership",
             "repTax", "write the taxonomy report for a genome source",
             "neighbors", "list the N closest neighboring representative genomes to all genomes in a source",
             "distDist", "find a set of genomes that are various distributed distances from a key genome",
             "buildFromGtos", "build a representative genome database for all genomes in a source",
             "hiCheck", "check the master representation list to see which genomes violate repgen hierarchy",
             "repSub", "validate the repgen sets to insure they are hierarchical",
             "buildSets", "create neighbor sets from a genome representation list file"
    };

    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        switch (command) {
        case "buildFromGtos" -> processor = new GtoRepBuildProcessor();
        case "gtoReps" -> processor = new GtoRepGenomeProcessor();
        case "ssuReps" -> processor = new SsuRepGenomeProcessor();
        case "distances" -> processor = new DistanceMatrixProcessor();
        case "repdb" -> processor = new RepGenomeDbProcessor();
        case "group" -> processor = new RepMatrixProcessor();
        case "classify" -> processor = new ClassifyProcessor();
        case "genomes" -> processor = new GenomeProcessor();
        case "md5" -> processor = new MD5Processor();
        case "coupling" -> processor = new CouplesProcessor();
        case "prepare" -> processor = new PrepareProcessor();
        case "seqTable" -> processor = new SeqTableProcessor();
        case "seqComp" -> processor = new SeqCompProcessor();
        case "seqTest" -> processor = new SeqTestProcessor();
        case "target" -> processor = new TargetProcessor();
        case "univ" -> processor = new UniRoleProcessor();
        case "missing" -> processor = new MissingRoleProcessor();
        case "gtoClass" -> processor = new GtoClassProcessor();
        case "ssuCheck" -> processor = new BadRnaProcessor();
        case "prio" -> processor = new PrioritizeProcessor();
        case "list" -> processor = new ListProcessor();
        case "rnaVerify" -> processor = new RnaVerifyProcessor();
        case "seedTable" -> processor = new SeedTableProcessor();
        case "build" -> processor = new BuildRepDbProcessor();
        case "buildRefDb" -> processor = new BuildRefDbProcessor();
        case "updateMaster" -> processor = new UpdateMasterProcessor();
        case "taxReport" -> processor = new TaxReportProcessor();
        case "repTax" -> processor = new RepTaxProcessor();
        case "neighbors" -> processor = new NeighborProcessor();
        case "distDist" -> processor = new DistributedDistanceProcessor();
        case "hiCheck" -> processor = new HierarchyCheckProcessor();
        case "repSub" -> processor = new RepGenSubCheckProcessor();
        case "buildSets" -> processor = new BuildSetsProcessor();
        case "-h", "--help" -> processor = null;
        default -> throw new RuntimeException("Invalid command " + command + ".");
        }
        if (processor == null)
            BaseProcessor.showCommands(COMMANDS);
        else {
            processor.parseCommand(newArgs);
            processor.run();
        }
    }
}
