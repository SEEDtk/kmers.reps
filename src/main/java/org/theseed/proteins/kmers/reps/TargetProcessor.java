/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.GenomeKmers;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command is designed to find a DNA target for a small set of genomes within the larger set.  The goal is to
 * find a DNA kmer that is in the seed protein gene for all of the small-set genomes but not in any of the large-set
 * genomes.  The positional parameters are the genome source (defaults to a GTO directory, but may be a file of PATRIC
 * IDs or a master directory) and GTO files for the small-set genomes.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -K	kmer size (default 20)
 *
 * --type	type of genome source (MASTER, DIR, or PATRIC)
 *
 * @author Bruce Parrello
 *
 */
public class TargetProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TargetProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** set of target genome IDs */
    private Set<String> targets;

    // COMMAND-LINE OPTIONS

    /** type of genome source (master directory, normal, PATRIC input) */
    @Option(name = "--type", usage = "type of input (master directory, GTO directory, PATRIC ID file)")
    private GenomeSource.Type inType;

    /** DNA kmer length */
    @Option(name = "-K", aliases = { "--kmer" }, metaVar = "21", usage = "length of kmer to target")
    private int kmerSize;

    /** genome source to process */
    @Argument(index = 0, metaVar = "gtoDir", usage = "genome directory to scan", required = true)
    private File gtoDir;

    /** small-set (target) genomes */
    @Argument(index = 1, metaVar = "gto1 gto2 ...", usage = "target genomes in set", required = true)
    private List<File> targetFiles;

    @Override
    protected void setDefaults() {
        this.inType = GenomeSource.Type.DIR;
        this.kmerSize = 20;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (this.kmerSize < 3)
            throw new ParseFailureException("Kmer size must be >= 3.");
        if (! this.gtoDir.exists())
            throw new FileNotFoundException("Genome source " + this.gtoDir + " is not found.");
        // Save the kmer size.
        DnaKmers.setKmerSize(this.kmerSize);
        GenomeKmers.setKmerSize(this.kmerSize);
        log.info("DNA kmer size is {}.", DnaKmers.kmerSize());
        // Get the set of target genomes.
        for (File targetFile : targetFiles) {
            if (! targetFile.canRead())
                throw new FileNotFoundException("Genome file " + targetFile + " is not found or unreadable.");
        }
        log.info("{} target genomes specified.", this.targetFiles.size());
        // Set up the genome source.
        this.genomes = this.inType.create(this.gtoDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.gtoDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // We must run through the target genomes to create a set of common DNA kmers.  We will collect the genome IDs
        // so we can skip them if the also appear in the source set.
        this.targets = new TreeSet<String>();
        Iterator<File> targetIter = this.targetFiles.iterator();
        // Start with the first target genome.
        File targetFile = targetIter.next();
        DnaKmers kmerSet = this.getSeedKmers(targetFile);
        // Take the intersection with the other targets.
        while (targetIter.hasNext()) {
            DnaKmers otherSet = this.getSeedKmers(targetIter.next());
            kmerSet.retainAll(otherSet);
        }
        log.info("{} kmers found in target seed proteins.", kmerSet.size());
        // Now run through the genomes from the source.
        int gCount = 0;
        Iterator<Genome> iter = this.genomes.iterator();
        while (kmerSet.size() > 0 && iter.hasNext()) {
            Genome genome = iter.next();
            if (this.targets.contains(genome.getId()))
                log.info("Skipping target genome {}.", genome);
            else {
                gCount++;
                log.info("Processing genome {}: {}.", gCount, genome);
                GenomeKmers gKmers = new GenomeKmers(genome);
                kmerSet.removeAll(gKmers);
                log.info("{} kmers after processing genome {}.", kmerSet.size(), genome);
            }
        }
        // Write the kmers found to the output.
        if (kmerSet.size() <= 0)
            log.error("Could not find any eligible kmers.");
        else {
            log.info("Writing output.");
            for (String kmer : kmerSet)
                System.out.println(kmer);
        }
    }

    /**
     * Find the DNA kmer set for the seed protein in the specified file's genome.
     *
     * @param targetFile	file containing the genome
     *
     * @return the DNA kmers for the genome's seed protein
     *
     * @throws IOException
     */
    private DnaKmers getSeedKmers(File targetFile) throws IOException {
        Genome genome = new Genome(targetFile);
        log.info("Processing {} from {}.", genome, targetFile);
        this.targets.add(genome.getId());
        String fid = SeqTableProcessor.findSeed(genome);
        if (fid == null)
            throw new IllegalArgumentException("Genone " + genome.toString() + " has no valid seed protein.");
        DnaKmers retVal = new DnaKmers(genome.getDna(fid));
        return retVal;
    }

}
