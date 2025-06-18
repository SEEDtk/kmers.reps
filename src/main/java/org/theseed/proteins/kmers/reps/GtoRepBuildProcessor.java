package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.BaseGenomeProcessor;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.reps.RepGenomeDb.Representation;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * This class reads all the genomes from a genome source and builds a representative-genome database. It uses
 * the "stingy addition" method. That is, it adds a genome to the database only if it is outside the similarity
 * distance from all the genomes currently in the database.
 * 
 * The positional parameters are the genome source file or directory, the name of the role definition file, 
 * and the name of the output rep-genome database.  The role file should contain only a single role ID, but 
 * it may contain multiple definitions if the role has more than one name.
 *
 * The command-line option are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -t	genome source type (default DIR)
 * -K	kmer size (default 8)
 * -m	minimum similarity threshold for representation (default 100)
 * 
 * --verify     if specified, a FASTA file containing seed protein sequences for multiple genomes; a report
 *              will be produced of all the genomes that are not represented in the result database; the
 *              report will be written to the standard output
 * 
 * @author Bruce Parrello
 *
*/
public class GtoRepBuildProcessor extends BaseGenomeProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(GtoRepBuildProcessor.class);
    /** output representative-genome database */
    private RepGenomeDb repDb;
    /** seed role definitions */
    private RoleMap roles;
    /** number of rep-genomes created */
    private int repGenomeCount;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name = "--kmer", aliases = { "-K" }, metaVar = "9", usage = "kmer size")
    private int kmerSize;

    /** minimum similarity threshold */
    @Option(name = "--sim", aliases = { "-m" }, metaVar = "150", usage = "minimum similarity threshold")
    private int minSim;

    /** verification FASTA file */
    @Option(name = "--verify", metaVar = "seeds.fasta",
            usage = "if specified, a FASTA file containing seed protein sequences for multiple genomes; a report will be produced of all the genomes that are not represented in the result database")
    private File verifyFile;

    /** seed role definition file */
    @Argument(index = 1, metaVar = "role.definitions", usage = "role definition file (should be 1 role only)", required = true)
    private File roleFile;

    /** output database name */
    @Argument(index = 2, metaVar = "repDb", usage = "output representative-genome database name", required = true)
    private File outFile;

    @Override
    protected void setSourceDefaults() {
        this.setLevel(Details.PROTEINS);
        this.kmerSize = 8;
        this.minSim = 100;
        this.verifyFile = null;
    }

    @Override
    protected void validateSourceParms() throws IOException, ParseFailureException {
        if (! this.roleFile.canRead())
            throw new IOException("Role definition file is not found or unreadable.");
        if (this.verifyFile != null && ! this.verifyFile.canRead())
            throw new IOException("Verification file is not found or unreadable.");
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer size must be at least 2.");
        if (this.minSim < 1)
            throw new ParseFailureException("Minimum similarity must be at least 1.");
        // Save the kmer size.
        ProteinKmers.setKmerSize(this.kmerSize);
        log.info("Kmer size is {}. Minimum similarity is {}.", this.kmerSize, this.minSim);
        // Read the role definitions.
        this.roles = RoleMap.load(this.roleFile);
        if (this.roles.size() != 1)
            throw new ParseFailureException("Role definition file must contain exactly 1 role.");
        // Insure we can write to the output file.
        try (PrintWriter writer = new PrintWriter(this.outFile)) {
            writer.println("");
        }
        // Get an array of the role names.
        String[] roleNames = this.roles.values().toArray(String[]::new);
        // Initialize the representative-genome database.
        this.repDb = new RepGenomeDb(this.minSim, roleNames);
    }

    @Override
    protected void runCommand() throws Exception {
        // We loop through the genome source, creating a RepGenome for each incoming genome, and
        // then check it against the current repDb.
        Collection<String> genomeIds = this.getGenomeIds();
        int gTotal = genomeIds.size();
        log.info("Scanning {} genomes from source.", gTotal);
        long lastMsg = System.currentTimeMillis();
        List<RepGenome> repGenomes = genomeIds.parallelStream().map(x -> this.computeRepGenome(x)).collect(Collectors.toList());
        int gCount = 0;
        gTotal = repGenomes.size();
        for (RepGenome repGenome : repGenomes) {
            gCount++;
            this.repDb.checkGenome(repGenome);
            if (System.currentTimeMillis() - lastMsg >= 10000) {
                lastMsg = System.currentTimeMillis();
                log.info("{} of {} genomes processed, {} in database.", gCount, gTotal, this.repDb.size());
            }
        }
        log.info("{} genomes in the representative-genome database.", this.repDb.size());
        // Save the representative-genome database.
        this.repDb.save(this.outFile);
        log.info("Wrote rep-genome database to {}.", this.outFile);
        // Verify the result if requested.
        if (this.verifyFile == null)
            System.out.println("No verification requested.");
        else {
            // Here we would verify the representative-genome database against the verification file.
            System.out.println("seq_id\tname\tsimilarity");
            log.info("Verifying repDb in {} against {}.", this.outFile, this.verifyFile);
            lastMsg = System.currentTimeMillis();
            int badCount = 0;
            int count = 0;
            try (FastaInputStream fastaStream = new FastaInputStream(this.verifyFile)) {
                for (Sequence seq : fastaStream) {
                    count++;
                    // Find the closest representative for this sequence.
                    Representation rep = this.repDb.findClosest(seq);
                    if (! rep.isRepresented()) {
                        // Here we have to output the bad genome.
                        System.out.format("%s\t%s\t%d%n", seq.getLabel(), seq.getComment(), rep.getSimilarity());
                        badCount++;
                    }
                    if (System.currentTimeMillis() - lastMsg >= 10000) {
                        lastMsg = System.currentTimeMillis();
                        log.info("{} sequences processed. {} bad found.", count, badCount);
                    }
                }
            }
        }
    }

    /**
     * This method searches for the seed protein in the genome and creates a representative-genome summary
     * for it. If multiple copies of the seed protein are found, the longest one is kept.
     * 
     * @param genomeId    ID of the genome to search
     */
    private RepGenome computeRepGenome(String genomeId) {
        String bestProtein = "";
        int protLen = 0;
        String fid = null;
        Genome genome = this.getGenome(genomeId);
        for (Feature feat : genome.getFeatures()) {
            if (feat.isProtein() && feat.isInteresting(this.roles)) {
                // Here we have a seed protein.
                String prot = feat.getProteinTranslation();
                // Remember the longest protein string.
                if (prot.length() > protLen) {
                    bestProtein = prot;
                    protLen = prot.length();
                    fid = feat.getId();
                }
            }
        }
        RepGenome retVal;
        if (fid == null) {
            // Here no seed protein was found.
            log.warn("No seed protein found for {}.", genome);
            retVal = null;
        } else {
            // Here we can return a genome representation object.
            retVal = new RepGenome(fid, genome.getName(), bestProtein);
        }
        synchronized (this) {
            if (retVal != null) {
                this.repGenomeCount++;
                if (this.repGenomeCount % 100 == 0)
                    log.info("Computed {} representative genomes.", this.repGenomeCount);
            }
        }
        return retVal;
    }

}
