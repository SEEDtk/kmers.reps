/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.GenomeDescriptor;
import org.theseed.sequence.GenomeDescriptorSet;
import org.theseed.sequence.GenomeDescriptorSet.Rating;
import org.theseed.sequence.GenomeDescriptorSource;
import org.theseed.sequence.ProteinKmers;

/**
 * This command compares a set of genomes in a genome source to a four-column table produced by the "seqTable" command-- the
 * "reference set".  For each genome from the source, the command finds the closest genome in the reference set using
 * various criteria.  A report is produced allowing comparison of the results.
 *
 * This program can also optionally produce an outlier report.  In the outlier report, each mismatch
 * is shown along with its protein similarity score (seed_seed_sim) and SSU similarity score
 * (rna_rna_sim).  We then display metrics for various cross-matches. The "rna_seed_sim" indicates
 * how close the SSU-closest genome is to the test genome using the protein similarity score.  If
 * this is above the RepGen threshold, then the test genome is in range of both representatives,
 * and the mismatch is minor.  "reps_seed_sim" indicates the protein similarity between the
 * two representatives.  If this is close to the threshold, it is another indicator the mismatch
 * is not too bad.  Similar information is provided using SSU matches.  "seed_rna_sim" is the
 * SSU similarity between the test genome and the seed representative.  "reps_rna_sim" is the
 * SSU similarity between the seed representative and the SSU representative.  These latter
 * numbers will be more useful once we find the SSU threshold that most closely matches the protein
 * threshold.
 *
 * The positional parameters are the name of the reference set file and the name of the testing genome source.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output report file (if not STDOUT)
 * -t	type of genome source (default DIR)
 * -x	if specified, the testing genome source is a four-column table instead of a genome source
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
    private static final Logger log = LoggerFactory.getLogger(SeqTestProcessor.class);
    /** reference genome set */
    private GenomeDescriptorSet refGenomes;
    /** genome descriptor resource, if any */
    private AutoCloseable descriptorResource;
    /** testing genome source */
    private Iterator<GenomeDescriptor> testGenomes;
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

    /** if specified, the input is presumed to be a four-column table instead of a genome source */
    @Option(name = "--4col", aliases = { "-4", "-x", "--seqs" }, usage = "source is a four-column table file")
    private boolean fourColFlag;

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
        this.fourColFlag = false;
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
        if (this.fourColFlag) {
            // Here we have a four-column table.
            log.info("Reading from four-column table at {}.", this.testGenomesDir);
            GenomeDescriptor.FileIter descIter = new GenomeDescriptor.FileIter(this.testGenomesDir);
            this.descriptorResource = descIter;
            this.testGenomes = descIter;
        } else {
            // Here we have a standard genome source.
            log.info("Data will be extracted from genome source {}.", this.testGenomesDir);
            GenomeDescriptorSource source = new GenomeDescriptorSource(this.testGenomesDir, sourceType);
            this.testGenomes = source.iterator();
        }
        // Insure the reference-genome file exists.
        if (! this.refGenomesFile.canRead())
            throw new IOException("Reference-genome file " + this.refGenomesFile + " is not found or unreadable.");
        // Open the outlier report if it is specified.
        if (this.outlierFile != null) {
            this.outlierWriter = new PrintWriter(this.outlierFile);
            log.info("Outlier report will be written to {}.", this.outlierFile);
            // Write the report header.
            this.outlierWriter.println("genome_id\tname\tseed_ref_id\tseed_ref_name\tseed_seed_sim\tseed_rna_sim\trna_ref_id\trna_ref_name\trna_rna_sim\trna_seed_sim\treps_seed_sim\treps_rna_sim");
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
            writer.println("genome_id\tname\tseed_rep_id\tseed_rep_sim\trna_rep_id\trna_rep_sim");
            // Here is where we will store our results.
            GenomeDescriptorSet.Rating seedResult;
            GenomeDescriptorSet.Rating rnaResult;
            // Now we loop through the testing genomes.
            int count = 0;
            int goodSims = 0;
            long start = System.currentTimeMillis();
            while (this.testGenomes.hasNext()) {
                GenomeDescriptor testDescriptor = this.testGenomes.next();
                count++;
                log.debug("Processing genome {}: {}.", count, testDescriptor);
                // Find the closest in each mode.
                seedResult = this.refGenomes.findClosest(testDescriptor,
                        GenomeDescriptorSet.FinderType.SEED_SIMILARITY);
                rnaResult =  this.refGenomes.findClosest(testDescriptor,
                        GenomeDescriptorSet.FinderType.RNA_SIMILARITY);
                if (seedResult.isSameGenome(rnaResult))
                    goodSims++;
                else {
                    log.debug("Sims mismatch for {}.", testDescriptor);
                    if (this.outlierWriter != null)
                        this.produceOutlierReport(testDescriptor, seedResult, rnaResult);
                }
                // Output the results.  Each genome is an output line, with the results in sequence.
                writer.println(testDescriptor.getId() + "\t" + testDescriptor.getName() + "\t" +
                        seedResult.output() + "\t" + rnaResult.output());
                if (log.isInfoEnabled() && count % 100 == 0) {
                    double pct = goodSims * 100.0 / count;
                    double rate = count * 1000.0 / (System.currentTimeMillis() - start);
                    log.info(String.format("%d genomes processed.  %4.2f%% good, %4.2f genomes/second.",
                            count, pct, rate));
                }
            }
            log.info("{} of {} genomes got the same similarity results.", goodSims, count);
        } finally {
            if (this.outlierWriter != null)
                this.outlierWriter.close();
            if (this.descriptorResource != null)
                this.descriptorResource.close();
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
