/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.reports.CouplingReporter;
import org.theseed.utils.BaseProcessor;

/**
 * This command computes functional couplings using genomes in a specific genome directory.  It allows specification of
 * the method for computing the basis of the coupling and the type of reporting.  The default is to couple by
 * protein families and output a list of the genomes where the coupling is found.
 *
 * The positional parameter is the name of the directory containing the input genomes.
 *
 * The command-line options are as follows.
 *
 * -h	display command usage
 * -v	display more detailed status messages
 * -d	maximum acceptable distance for features to be considered neighbors (default 5000)
 * -t	classification method for features (default PROTFAM)
 * -m	minimum number of genomes for a pair to be output (default 10)
 *
 * --format		report format (default GROUP)
 *
 * @author Bruce Parrello
 *
 */
public class CouplesProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CouplesProcessor.class);
    /** hash of pairs to genome sets */
    private Map<FeatureClass.Pair, Set<String>> pairMap;
    /** feature classifier */
    private FeatureClass classifier;

    // COMMAND-LINE OPTIONS

    /** maximum gap distance */
    @Option(name = "-d", aliases = { "--distance", "--gap" }, metaVar = "4500", usage = "maximum acceptable gap between coupled features")
    private int maxGap;

    /** classification method */
    @Option(name = "-t", aliases = { "--type" }, usage = "type of classification for features")
    private FeatureClass.Type classType;

    /** report type */
    @Option(name = "--format", aliases = { "--output", "--outFmt" }, usage = "output format")
    private CouplingReporter.Type reportType;

    /** minimum group size */
    @Option(name = "-m", aliases = { "--min" }, usage = "minimum group size for output to the report")
    private int minGroup;

    /** input directory */
    @Argument(index = 0, metaVar = "genomeDir", usage = "input genome directory", required = true)
    private File genomeDir;


    @Override
    protected void setDefaults() {
        this.classType = FeatureClass.Type.PGFAMS;
        this.reportType = CouplingReporter.Type.GROUP;
        this.minGroup = 10;
        this.maxGap = 5000;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (this.minGroup < 1)
            throw new IllegalArgumentException("Invalid minimum group size.  Must be at least 1.");
        if (this.maxGap < 1)
            throw new IllegalArgumentException("Invalid maximum gap size.  Must be at least 1.");
        if (! this.genomeDir.isDirectory())
            throw new FileNotFoundException("Specified genome directory" + this.genomeDir + " is not found or invalid.");
        // Create the feature classifier.
        this.classifier = this.classType.create();
        // Create the pair map.
        this.pairMap = new HashMap<FeatureClass.Pair, Set<String>>(100000);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        log.info("Initializing report.");
        // Start by creating the reporter object.
        try (CouplingReporter reporter = this.reportType.create(System.out, this.classifier)) {
            // Initialize the report.
            reporter.writeHeader();
            // Open the genome directory.
            log.info("Scanning genome directory {}.", this.genomeDir);
            GenomeDirectory genomes = new GenomeDirectory(this.genomeDir);
            int total = genomes.size();
            log.info("{} genomes found.", total);
            // Loop through the genomes, filling the pair map.
            int count = 0;
            for (Genome genome : genomes) {
                count++;
                log.info("Processing genome {} of {}: {}.", count, total, genome);
                List<FeatureClass.Result> gResults = this.classifier.getResults(genome);
                log.info("{} classifiable features found.", gResults.size());
                // Loop through the results.  For each one, we check subsequent results up to the gap
                // distance.
                int n = gResults.size();
                int pairCount = 0;
                for (int i = 0; i < n; i++) {
                    FeatureClass.Result resI = gResults.get(i);
                    for (int j = i + 1; j < n && resI.getDistance(gResults.get(j)) <= this.maxGap; j++) {
                        // Here we have two neighboring features, so we pair the classes.
                        for (String classI : resI) {
                            for (String classJ : gResults.get(j)) {
                                if (! classI.contentEquals(classJ)) {
                                    FeatureClass.Pair pair = this.classifier.new Pair(classI, classJ);
                                    Set<String> pairList = this.pairMap.computeIfAbsent(pair, k -> new HashSet<String>(5));
                                    pairList.add(genome.getId());
                                    pairCount++;
                                }
                            }
                        }
                    }
                }
                log.info("{} pairs found in {}.", pairCount, genome);
                // Register the genome with the report facility.
                reporter.register(genome);
            }
            // Now we produce the output.
            log.info("{} distinct pairs found in genome set.", this.pairMap.size());
            int outputCount = 0;
            int groupTotal = 0;
            int groupMax = 0;
            int skipCount = 0;
            count = 0;
            for (Map.Entry<FeatureClass.Pair, Set<String>> pairData : this.pairMap.entrySet()) {
                Set<String> gSet = pairData.getValue();
                // Only use the group if it is big enough.
                if (gSet.size() < this.minGroup) {
                    skipCount++;
                } else {
                    reporter.writePair(pairData.getKey(), gSet);
                    outputCount++;
                    if (gSet.size() > groupMax) groupMax = gSet.size();
                    groupTotal += gSet.size();
                }
                count++;
                if (log.isDebugEnabled() && count % 5000 == 0)
                    log.debug("{} pairs processed for output. {} skipped.", outputCount, skipCount);
            }
            // Finish the report.
            log.info("Finishing report. {} pairs output, {} skippped.", outputCount, skipCount);
            log.info("Mean group size is {}. Max group size is {}.", (double) groupTotal / outputCount, groupMax);
            reporter.finish();
        }
    }

}
