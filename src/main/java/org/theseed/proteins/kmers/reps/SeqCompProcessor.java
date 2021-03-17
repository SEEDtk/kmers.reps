/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This script is used to determine the correspondence between PheS distance and SSU-rRNA distance.  It takes an
 * output file from the SeqTable command and produces a two-column output file, each record containing the PheS
 * distance and SSU-rRNA distance for a single pair of genomes.  The log will display the Pearson's correlation
 * and other statistical measures.
 *
 * The input file is taken from the standard input.  There are no positional parameters.  The command-line options
 * are as follows.
 *
 * -h	show command-line usage
 * -v	show more detailed log messages
 * -K	protein kmer size (default 8)
 * -k	dna kmer size (default 15)
 * -i	name of the input file (if not STDIN)
 *
 * @author Bruce Parrello
 *
 */
public class SeqCompProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SeqCompProcessor.class);
    /** input stream */
    private TabbedLineReader inStream;

    // COMMAND-LINE OPTIONS

    /** input file name (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, usage = "name of the input file (if not STDIN)")
    private File inFile;

    /** protein kmer size */
    @Option(name = "-K", aliases = { "--protKmer" }, metaVar = "9", usage = "protein kmer size")
    private int protK;

    /** DNA kmer size */
    @Option(name = "-k", aliases = { "--dnaKmer" }, metaVar = "15", usage = "DNA kmer size")
    private int dnaK;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.protK = 8;
        this.dnaK = 15;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the kmers are reasonable.
        if (this.protK < 3)
            throw new ParseFailureException("Protein kmer size must be 3 or more.");
        if (this.dnaK < 4)
            throw new ParseFailureException("DNA kmer size must be 4 or more.");
        if (this.inFile == null) {
            log.info("Input will be from the standard input.");
            this.inStream = new TabbedLineReader(System.in);
        } else {
            log.info("Input will be from {}.", this.inFile);
            this.inStream = new TabbedLineReader(this.inFile);
        }
        // Set up the kmer sizes.
        DnaKmers.setKmerSize(this.dnaK);
        ProteinKmers.setKmerSize(this.protK);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Write the output header.
        System.out.println("prot_dist\tssu_dist");
        // We fill two lists with the input data.
        List<ProteinKmers> proteins = new ArrayList<ProteinKmers>(200);
        List<DnaKmers> ssuRnas = new ArrayList<DnaKmers>(200);
        log.info("Reading sequences from input file.");
        for (TabbedLineReader.Line line : this.inStream) {
            proteins.add(new ProteinKmers(line.get(1)));
            ssuRnas.add(new DnaKmers(line.get(2)));
        }
        this.inStream.close();
        // Set up to create the pairwise distances.
        int gCount = proteins.size();
        int pairs = gCount * (gCount + 1) / 2;
        log.info("{} genomes to process. {} pairs will be output.", gCount, pairs);
        double[] protDists = new double[pairs];
        double[] ssuDists = new double[pairs];
        log.info("Performing pairwise comparisons.");
        int pairCount = 0;
        // Loop through the possible pairs.
        long start = System.currentTimeMillis();
        for (int i = 0; i < gCount; i++) {
            for (int j = i + 1; j < gCount; j++) {
                protDists[pairCount] = proteins.get(i).distance(proteins.get(j));
                ssuDists[pairCount] = ssuRnas.get(i).distance(ssuRnas.get(j));
                System.out.format("%8.4f\t%8.4f%n", protDists[pairCount], ssuDists[pairCount]);
                pairCount++;
                if (log.isInfoEnabled() && (pairCount % 1000) == 0) {
                    Duration duration = Duration.ofMillis(System.currentTimeMillis() - start).dividedBy(pairCount);
                    log.info("{} of {} pairs computed, {} per pair.", pairCount, pairs, duration);
                }
            }
        }
        // All done.  Make sure the output is safe.
        System.out.flush();
        // Display the statistics.
        PearsonsCorrelation pComputer = new PearsonsCorrelation();
        double pCorr = pComputer.correlation(protDists, ssuDists);
        KendallsCorrelation kComputer = new KendallsCorrelation();
        double kCorr = kComputer.correlation(protDists, ssuDists);
        log.info("Pearson correlation = {}.  Kendall tau = {}.", pCorr, kCorr);
    }

}
