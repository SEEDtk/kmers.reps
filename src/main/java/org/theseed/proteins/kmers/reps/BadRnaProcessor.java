/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.CloseSets;
import org.theseed.reports.NeighborReporter;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.GenomeDescriptor;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseReportProcessor;

/**
 * This command attempts to identify bad SSU rRNA sequences. The basic technique is to find the closest genomes in a
 * reference set using both PheS and 16s RNA distances.  The closest genomes for each type of distance are then compared.
 *
 * The positional parameters are the name of the sequence file for the testing set and the name of the testing file for
 * the reference set.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output report file (if not STDOUT)
 *
 * --format		output report format
 * --dnaK		DNA kmer size (default 15)
 * --protK		protein kmer size (default 8)
 *
 * @author Bruce Parrello
 *
 */
public class BadRnaProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BadRnaProcessor.class);
    /** map of reference genomes to genome descriptors */
    private Map<String, GenomeDescriptor> referenceMap;
    /** output reporting object */
    private NeighborReporter reporter;
    /** estimated number of genomes for maps */
    private static final int ESTIMATED_GENOMES = 2000;


    // COMMAND-LINE OPTIONS

    /** output report format */
    @Option(name = "--format", usage = "type of report to produce")
    private NeighborReporter.Type reportType;

    /** DNA kmer size */
    @Option(name = "--dnaK", metaVar = "18", usage = "SSU rRNA kmer size")
    private int dnaK;

    /** protein kmer size */
    @Option(name = "--protK", metaVar = "9", usage = "seed protein kmer size")
    private int protK;

    /** testing set sequences */
    @Argument(index = 0, metaVar = "testing.tbl", usage = "file of identifying sequences for genomes to test", required = true)
    private File testingFile;

    /** reference set sequences */
    @Argument(index = 1, metaVar = "reference.tbl", usage = "file of identifying sequences for reference genomes", required = true)
    private File referenceFile;

    @Override
    protected void setReporterDefaults() {
        this.dnaK = 15;
        this.protK = 8;
        this.reportType = NeighborReporter.Type.STATS;
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
        // Verify the testing set.
        if (! this.testingFile.canRead())
            throw new FileNotFoundException("Testing set file " + this.testingFile + " is not found or invalid.");
        // Read in the reference set.
        this.referenceMap = new HashMap<String, GenomeDescriptor>(ESTIMATED_GENOMES);
        log.info("Reading reference set from {}.", this.referenceFile);
        try (GenomeDescriptor.FileIter iter = new GenomeDescriptor.FileIter(this.referenceFile)) {
            while (iter.hasNext()) {
                GenomeDescriptor newGenome = iter.next();
                this.referenceMap.put(newGenome.getId(), newGenome);
            }
        }
        log.info("{} genomes found in reference set.", this.referenceMap.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        this.reporter = this.reportType.create(writer);
        this.reporter.openReport();
        int count = 0;
        // Loop through the testing file, forming neighbor lists for each genome.
        try (GenomeDescriptor.FileIter iter = new GenomeDescriptor.FileIter(this.testingFile)) {
            while (iter.hasNext()) {
                // Get the current testing genome and create a result repository for it.
                count++;
                GenomeDescriptor current = iter.next();
                CloseSets results = new CloseSets();
                log.info("Processing genome #{}: {} ({}).", count, current.getId(), current.getName());
                // Compare the genome to everything in the reference set.
                for (GenomeDescriptor refGenome : this.referenceMap.values()) {
                    double protDist = current.getSeedDistance(refGenome);
                    double rnaDist = current.getRnaDistance(refGenome);
                    results.addGenome(refGenome.getId(), protDist, rnaDist);
                }
                // Output the results.
                this.reporter.recordGenome(current.getId(), current.getName(), results);
            }
        }
        // Finish the report.
        this.reporter.finish();
    }

}
