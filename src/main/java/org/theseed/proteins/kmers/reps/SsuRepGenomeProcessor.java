/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.SsuGenome;
import org.theseed.sequence.SsuRepresentation;

/**
 * This command computes the closest representative genome to each input genome using SSUs from a
 * seedTable file.  In this file, the representative genome ID is in a column called "genome_id",
 * its name is in a column called "genome_name", and its SSU RNA sequence is in a column called
 * "ssu_rna".  The SSU from each input genome will be extracted and compared to the SSUs from the
 * seed-table file to determine the closest.  Because SSUs vary widely in size, distances will be
 * used instead of similarity counts.
 *
 * The positional parameters are the name of the seed-table file and the name of the input genome
 * source for the genomes to process.  The output report will be to the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 * -K	kmer size to use (default 18)
 * -t	type of genome source (default DIR)
 *
 * @author Bruce Parrello
 *
 */
public class SsuRepGenomeProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SsuRepGenomeProcessor.class);
    /** list of SSU genome objects to use to find the closest genome */
    private Set<SsuGenome> ssuList;
    /** genome source for input genomes */
    private GenomeSource genomes;

    // COMMAND-LINE OPTIONS

    /** kmer size to use */
    @Option(name = "kmer", aliases = { "-K", "--kmerSize" }, metaVar = "20", usage = "DNA kmer size to use")
    private int kmerSize;

    /** type of genome source */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of input genome source")
    private GenomeSource.Type sourceType;

    /** input seed-table file */
    @Argument(index = 0, metaVar = "ssuTable.tbl", usage = "file containing genome IDs, names, and SSUs",
            required = true)
    private File ssuFile;

    /** input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome directory or file")
    private File inDir;

    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.kmerSize = 18;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify and set the kmer size.
        if (this.kmerSize < 3)
            throw new ParseFailureException("Kmer size must be at least 3.");
        DnaKmers.setKmerSize(this.kmerSize);
        log.info("Kmer size is {}.", this.kmerSize);
        // Verify the SSU file.
        if (! this.ssuFile.canRead())
            throw new FileNotFoundException("SSU input file " + this.ssuFile + " is not found or unreadable.");
        // Verify the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " does not exist.");
        log.info("Connecting to genome source type {} in {}.", this.sourceType, this.inDir);
        this.genomes = this.sourceType.create(this.inDir);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Load the SSUs.
        this.ssuList = new HashSet<>(2000);
        log.info("Loading SSUs from {}.", this.ssuFile);
        try (TabbedLineReader ssuStream = new TabbedLineReader(this.ssuFile)) {
            int ssuCol = ssuStream.findField("ssu_rna");
            int nameCol = ssuStream.findField("genome_name");
            int idCol = ssuStream.findField("genome_id");
            int count = 0;
            for (TabbedLineReader.Line line : ssuStream) {
                SsuGenome ssuGenome = new SsuGenome(line.get(idCol), line.get(nameCol), line.get(ssuCol));
                this.ssuList.add(ssuGenome);
                count++;
                if (log.isInfoEnabled() && count % 1000 == 0)
                    log.info("{} sequences stored.", count);
            }
        }
        log.info("{} representatives found in SSU file.", this.ssuList.size());
        // Write the output header.
        writer.println("genome_id\tgenome_name\trep_id\trep_name\tsimilarity\tdistance");
        // Loop through the input genomes, processing them.
        int count = 0;
        int errCount = 0;
        int notFoundCount = 0;
        for (Genome genome : this.genomes) {
            SsuGenome inputSsu = SsuGenome.create(genome);
            count++;
            if (inputSsu == null) {
                log.warn("No SSU found in {}.", genome);
                errCount++;
            } else {
                log.info("Processing genome #{}: {}.", count, genome);
                // Compare the genome to each representative, remembering the best.
                SsuRepresentation rep = new SsuRepresentation();
                for (SsuGenome candidate : this.ssuList)
                    rep.checkSsu(inputSsu, candidate);
                if (rep.isEmpty()) {
                    log.warn("No representative found for {}.", genome);
                    notFoundCount++;
                } else {
                    // Write out the results.
                    writer.format("%s\t%s\t%s\t%s\t%d\t%6.4f%n", inputSsu.getGenomeId(), inputSsu.getName(),
                            rep.getGenomeId(), rep.getGenomeName(), rep.getHitCount(), rep.getDistance());
                }
            }
        }
        log.info("{} genomes processed.  {} missing SSUs. {} outliers.", count, errCount, notFoundCount);
    }

}
