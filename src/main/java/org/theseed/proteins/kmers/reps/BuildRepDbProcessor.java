/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command builds a representative-genome database with pre-selected representatives and a
 * specific universal protein for the seed.  If the protein does not exist in the representative
 * genome, it will be skipped.
 *
 * The positional parameters are the name of the genome source for the pre-selected representatives,
 * the name of the role definition file, and the name of the output rep-genome database.  The role file
 * should contain only a single role ID, but it may contain multiple defintions if the role has more
 * than one name.
 *
 * The command-line options are as follows.
 *
 * -h	show command-line usage
 * -v	display more frequent log messages
 * -K	kmer size (default 8)
 * -t	genome source type (derfault DIR)
 *
 * @author Bruce Parrello
 *
 */
public class BuildRepDbProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BuildRepDbProcessor.class);
    /** rep-genome database being built */
    private RepGenomeDb repDb;
    /** source for input genomes */
    private GenomeSource genomes;
    /** list of role names */
    private String[] roleNames;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name = "--kmer", aliases = { "--kmerSize", "-K" }, metaVar = "9", usage = "protein kmer size")
    private int kmerSize;

    /** genome source type */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** role definitions */
    @Argument(index = 0, metaVar = "role.definitions", usage = "role definition file (should be 1 role only)",
            required = true)
    private File roleFile;

    /** genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "genome source directory or file", required = true)
    private File inDir;

    /** output file for repDB */
    @Argument(index = 2, metaVar = "outDb.ser", usage = "output file for rep-genome database", required = true)
    private File outFile;

    @Override
    protected void setDefaults() {
        this.kmerSize = 8;
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the kmer size.
        if (this.kmerSize < 3)
            throw new ParseFailureException("Kmer size must be at least 3.");
        ProteinKmers.setKmerSize(this.kmerSize);
        log.info("Kmer size is {}.", this.kmerSize);
        // Verify that the role definition file exists and get the role names.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or unreadable.");
        try (TabbedLineReader roleStream = new TabbedLineReader(this.roleFile, 3)) {
            List<String> names = new ArrayList<String>();
            var ids = new TreeSet<String>();
            for (TabbedLineReader.Line line : roleStream) {
                String id = line.get(0);
                ids.add(id);
                if (ids.size() > 1)
                    throw new IOException("More than one role ID found in " +  this.roleFile + ".");
                names.add(line.get(2));
            }
            if (names.size() < 1)
                throw new ParseFailureException("No role names found in " + this.roleFile + ".");
            this.roleNames = names.stream().toArray(String[]::new);
            log.info("{} names presented for role {}.", names.size(), ids.first());
        }
        // Load the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " not found.");
        this.genomes = this.sourceType.create(this.inDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.inDir);
        // Finally, verify that we can write to the output file.
        if (this.outFile.exists() && ! this.outFile.canWrite())
            throw new IOException("Cannot write to output file " + this.outFile + ".");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the rep-genome database.
        this.repDb = new RepGenomeDb(this.kmerSize, this.roleNames);
        // Loop through the genome source.  For each genome, add it to the repgenome database.
        int count = 0;
        final int nGenomes = this.genomes.size();
        int badGenomes = 0;
        int goodGenomes = 0;
        for (Genome genome : this.genomes) {
            count++;
            log.info("Processing genome {} of {}: {}.", count, nGenomes, genome);
            RepGenome rep = this.repDb.getSeedProtein(genome);
            if (rep == null) {
                log.warn("No seed protein found for {}.", genome);
                badGenomes++;
            } else {
                this.repDb.addRep(rep);
                goodGenomes++;
            }
        }
        // Now we save the results.
        log.info("Saving rep-genome database to {}.", this.outFile);
        this.repDb.save(this.outFile);
        log.info("{} genomes processed, {} stored, {} were missing seed proteins.", count, goodGenomes, badGenomes);
    }

}
