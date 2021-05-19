/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDescriptor;
import org.theseed.genome.GenomeDescriptorSet;
import org.theseed.genome.GenomeDescriptorSet.Rating;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command compares a set of genomes in a genome source to a four-column table produced by the "seqTable" command-- the
 * "reference set".  For each genome from the source, the command finds the closest genome in the reference set using
 * various criteria.  A report is produced allowing comparison of the results.
 *
 * The positional parameters are the name of the reference set file and the name of the testing genome source.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output report file (if not STDOUT)
 * -t	type of genome source (default DIR)
 *
 * --dnaK		DNA kmer size (default 15)
 * --protK		protein kmer size (default 8)
 * --outliers	name of a file to contain a report on mismatches
 *
 * @author Bruce Parrello
 *
 */
public class SeqTestProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SeqTestProcessor.class);
    /** reference genome set */
    private GenomeDescriptorSet refGenomes;
    /** testing genome source */
    private GenomeSource testGenomes;
    /** outlier report stream */
    private PrintWriter outlierWriter;

    // COMMAND-LINE OPTIONS

    /** genome source type */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** DNA kmer size */
    @Option(name = "--dnaK", metaVar = "18", usage = "SSU rRNA kmer size")
    private int dnaK;

    /** protein kmer size */
    @Option(name = "--protK", metaVar = "9", usage = "seed protein kmer size")
    private int protK;

    /** outliers report file */
    @Option(name = "--outliers", metaVar = "outliers.report.tbl", usage = "optional outliers report output file")
    private File outlierFile;

    /** reference set file */
    @Argument(index = 0, metaVar = "refGenomes.tbl", usage = "sequence table for reference genomes")
    private File refGenomesFile;

    /** testing set source */
    @Argument(index = 1, metaVar = "testGenomes", usage = "genome source (file or directory) for testing genomes")
    private File testGenomesDir;

    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.dnaK = 15;
        this.protK = ProteinKmers.kmerSize();
        this.outlierFile = null;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify the kmer sizes.
        if (this.dnaK < 3)
            throw new ParseFailureException("DNA kmer size must be >= 3.");
        if (this.protK < 1)
            throw new ParseFailureException("Protein kmer size must be >= 1.");
        // Store the kmer sizes.
        ProteinKmers.setKmerSize(this.protK);
        DnaKmers.setKmerSize(this.dnaK);
        log.info("DNA kmer size is {}.  Protein kmer size is {}.", this.dnaK, this.protK);
        // Set up the genome source.
        this.testGenomes = this.sourceType.create(this.testGenomesDir);
        log.info("{} genomes in testing set read from {}.", this.testGenomes.size(), this.testGenomesDir);
        // Insure the reference-genome file exists.
        if (! this.refGenomesFile.canRead())
            throw new IOException("Reference-genome file " + this.refGenomesFile + " is not found or unreadable.");
        // Open the outlier report if it is specified.
        if (this.outlierFile != null) {
            this.outlierWriter = new PrintWriter(this.outlierFile);
            log.info("Outlier report will be written to {}.", this.outlierFile);
            // Write the report header.
            this.outlierWriter.println("genome_id\tname\tseed_ref_id\tseed_ref_name\tseed_seed_sim\tseed_rna_sim\trna_ref_id\trna_ref_name\trna_seed_sim\trna_rna_sim\treps_seed_sim\treps_rna_sim");
        } else {
            this.outlierWriter = null;
            log.info("No outlier report will be written.");
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception  {
        try {
            // Read in the reference genomes.
            log.info("Reading reference genomes from {}.", this.refGenomesFile);
            this.refGenomes = new GenomeDescriptorSet(this.refGenomesFile);
            // Start the report.  For each finder type, we have a reference ID and a rating.
            writer.println("genome_id\tname\t" + Stream.of(GenomeDescriptorSet.FinderType.values()).map(x -> (x + "_ref_id\t" + x + "_rating"))
                    .collect(Collectors.joining("\t")));
            // Now we loop through the testing genomes.
            int count = 0;
            int goodSims = 0;
            int goodDist = 0;
            int goodSeed = 0;
            for (Genome testGenome : this.testGenomes) {
                count++;
                log.info("Processing genome {} of {}: {}.", count, this.testGenomes.size(), testGenome);
                // Get this genome's identifying sequences.
                GenomeDescriptor testDescriptor = new GenomeDescriptor(testGenome);
                // Find the closest in each mode.
                List<GenomeDescriptorSet.Rating> results = new ArrayList<>();
                for (GenomeDescriptorSet.FinderType type : GenomeDescriptorSet.FinderType.values()) {
                    GenomeDescriptorSet.Rating result = this.refGenomes.findClosest(testDescriptor, type);
                    results.add(result);
                }
                if (GenomeDescriptorSet.Rating.test(results, GenomeDescriptorSet.SIM_TYPES))
                    goodSims++;
                else {
                    log.warn("Sims mismatch for {}.", testGenome);
                    if (this.outlierWriter != null)
                        this.produceOutlierReport(testDescriptor,
                                results.get(GenomeDescriptorSet.FinderType.SEED_SIMILARITY.ordinal()),
                                results.get(GenomeDescriptorSet.FinderType.RNA_SIMILARITY.ordinal()));
                }
                if (GenomeDescriptorSet.Rating.test(results, GenomeDescriptorSet.DIST_TYPES))
                    goodDist++;
                else
                    log.warn("Distance mismatch for {}.", testGenome);
                if (GenomeDescriptorSet.Rating.test(results, GenomeDescriptorSet.SEED_TYPES))
                    goodSeed++;
                else
                    log.warn("Seed mismatch for {}.", testGenome);
                // Output the results.  Each genome is an output line, with the results in sequence.
                writer.println(testGenome.getId() + "\t" + testGenome.getName() + "\t" +
                        results.stream().map(x -> x.output()).collect(Collectors.joining("\t")));
            }
            log.info("{} of {} genomes got the same similarity results.", goodSims, count);
            log.info("{} of {} genomes got the same distance results.", goodDist, count);
            log.info("{} of {} genomes got the same seed-protein results.", goodSeed, count);
        } finally {
            if (this.outlierWriter != null)
                this.outlierWriter.close();
        }
    }

    /**
     * Process a similarity mismatch to produce an outlier report.
     *
     * @param testGenome	descriptor of the genome that produced the mismatch
     * @param resultSeed	closest-genome result using SEED protein
     * @param resultRna		closest-genome result using RNA sequence
     */
    private void produceOutlierReport(GenomeDescriptor testGenome, Rating resultSeed, Rating resultRna) {
        // Compute the converse proximities.  We know the seed-similarity to the seed-closest.  Find the RNA similarity of the seed-closest.
        double seedRnaSim = testGenome.getRnaSim(resultSeed.getGenome());
        // Similarly, we know the RNA similarity of the RNA-closest.  Find the seed similarity of the RNA-closest.
        double rnaSeedSim = testGenome.getSeedSim(resultRna.getGenome());
        // Get the proximities for the two reference genomes.
        double refSeedSim = resultSeed.getGenome().getSeedSim(resultRna.getGenome());
        double refRnaSim = resultSeed.getGenome().getRnaSim(resultRna.getGenome());
        // Write the output line.
        this.outlierWriter.format("%s\t%s\t%s\t%s\t%8.4f\t%8.4f\t%s\t%s\t%8.4f\t%8.4f\t%8.4f\t%8.4f%n",
                testGenome.getId(), testGenome.getName(), resultSeed.getGenomeId(), resultSeed.getGenomeName(),
                resultSeed.getProximity(), seedRnaSim, resultRna.getGenomeId(), resultRna.getGenomeName(),
                resultRna.getProximity(), rnaSeedSim, refSeedSim, refRnaSim);
    }

}
