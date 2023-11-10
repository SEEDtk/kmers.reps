/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.reports.ErrorSummaryReport;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.GenomeDescriptor;
import org.theseed.sequence.GenomeDescriptorSet;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseReportProcessor;

/**
 * This script is used to determine the correspondence between PheS distance and SSU-rRNA distance.  It takes an
 * output file from the SeqTable command and produces a two-column output file, each record containing the PheS
 * distance and SSU-rRNA distance for a single pair of genomes.  The log will display the Pearson's correlation
 * and other statistical measures.
 *
 * The positional parameter is the name of the input file.  The command-line options
 * are as follows.
 *
 * -h	show command-line usage
 * -v	show more detailed log messages
 * -K	protein kmer size (default 8)
 * -k	dna kmer size (default 15)
 * -o	output file (if not STDOUT)
 *
 * --summary	name of a file to contain a summary report that displays the mean and standard deviation of the
 * 				distance error for each genome
 * --summType	type of summary report
 *
 * @author Bruce Parrello
 *
 */
public class SeqCompProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SeqCompProcessor.class);
    /** genome descriptors from the input file */
    private GenomeDescriptorSet refGenomes;
    /** error summary reporting object */
    private ErrorSummaryReport summaryReporter;

    // COMMAND-LINE OPTIONS

    /** protein kmer size */
    @Option(name = "-K", aliases = { "--protKmer" }, metaVar = "9", usage = "protein kmer size")
    private int protK;

    /** DNA kmer size */
    @Option(name = "-k", aliases = { "--dnaKmer" }, metaVar = "15", usage = "DNA kmer size")
    private int dnaK;

    /** summary report file (optional) */
    @Option(name = "--summary", metaVar = "summaryReport.tbl", usage = "if specified, the name of a file to contain an error summary report")
    private File summaryFile;

    /** summary report type */
    @Option(name = "--summType", usage = "type of summary report (if a summary report is requested)")
    private ErrorSummaryReport.Type summType;

    /** input file name (if not STDIN) */
    @Argument(index = 0, metaVar = "rep.seqs.tbl", usage = "RepGen file containing PheS and SSU rRNA sequences", required = true)
    private File inFile;


    @Override
    protected void setReporterDefaults() {
        this.protK = 8;
        this.dnaK = 15;
        this.summaryFile = null;
        this.summType = ErrorSummaryReport.Type.STATS;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Insure the kmers are reasonable.
        if (this.protK < 3)
            throw new ParseFailureException("Protein kmer size must be 3 or more.");
        if (this.dnaK < 4)
            throw new ParseFailureException("DNA kmer size must be 4 or more.");
        // Insure the input file exists.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        // Read in the representative genomes.
        this.refGenomes = new GenomeDescriptorSet(this.inFile);
        // Set up the kmer sizes.
        DnaKmers.setKmerSize(this.dnaK);
        ProteinKmers.setKmerSize(this.protK);
        // Initialize the summary report (if any).
        this.summaryReporter = ErrorSummaryReport.create(this.summType, this.summaryFile);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the output header.
        writer.println("genome1\tname1\tgenome2\tname2\tseed_dist\trna_dist");
        // Compute the genome counts.
        int gCount = this.refGenomes.size();
        int pairs = gCount * (gCount + 1) / 2;
        log.info("{} genomes to process. {} pairs will be output.", gCount, pairs);
        double[] protDists = new double[pairs];
        double[] ssuDists = new double[pairs];
        log.info("Performing pairwise comparisons.");
        int pairCount = 0;
        // Loop through the genomes
        int count = 0;
        long start = System.currentTimeMillis();
        for (GenomeDescriptor genome1 : this.refGenomes) {
            count++;
            log.info("Processing genome {} of {}: {},", count, gCount, genome1);
            this.summaryReporter.registerGenome(genome1.getId(), genome1.getName());
            // Loop through all subsequent genomes.
            Iterator<GenomeDescriptor> iter = this.refGenomes.tailIter(genome1.getId());
            while (iter.hasNext()) {
                GenomeDescriptor genome2 = iter.next();
                protDists[pairCount] = genome1.getSeedDistance(genome2);
                ssuDists[pairCount] = genome1.getRnaDistance(genome2);
                this.summaryReporter.recordDistances(genome1.getId(), genome2.getId(), protDists[pairCount],
                        ssuDists[pairCount]);
                writer.format("%s\t%s\t%s\t%s\t%6.4f\t%6.4f%n", genome1.getId(), genome1.getName(), genome2.getId(),
                        genome2.getName(), protDists[pairCount], ssuDists[pairCount]);
                pairCount++;
                if (log.isInfoEnabled() && (pairCount % 1000) == 0) {
                    // The rate computation is complicated by the fact our times are in milliseconds and we want to report
                    // on pairs per second.  We multiply by 1000, but we must use a long integer to prevent overflow.
                    long rate = pairCount * 1000L / (System.currentTimeMillis() - start);
                    log.info("{} of {} pairs computed, {} pairs/second.", pairCount, pairs, rate);
                }
            }
        }
        // Display the summary report and the statistics.
        this.summaryReporter.writeReport();
        PearsonsCorrelation pComputer = new PearsonsCorrelation();
        double pCorr = pComputer.correlation(protDists, ssuDists);
        KendallsCorrelation kComputer = new KendallsCorrelation();
        double kCorr = kComputer.correlation(protDists, ssuDists);
        log.info("Pearson correlation = {}.  Kendall tau = {}.", pCorr, kCorr);
    }

}
