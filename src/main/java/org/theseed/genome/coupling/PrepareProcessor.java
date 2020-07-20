/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Coupling;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.TabbedLineReader;

/**
 * This command creates specially-augmented GTOs that can be used to drive the molecular machine processor.  Each
 * GTO will contain coupling information taken from the output of the CouplesProcessor with the SCORES option.
 *
 * The positional parameters are the name of the CouplesProcessor output file and the name of an input directory
 * containing GTOs.  The GTOs will be updated in place to contain coupling data.
 *
 * The command-line options are as follows:
 *
 * -h	show command-line usage
 * -v	display more detailed messages on the log
 * -t	type of feature classification used in the input
 * -n	algorithm for determining the feature neighborhood
 * -d	maximum acceptable distance for features to be considered neighbors
 *
 * @author Bruce Parrello
 *
 */
public class PrepareProcessor extends BaseCouplingProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PrepareProcessor.class);
    /** input directory */
    private GenomeDirectory inputGenomes;
    /** coupling map, keyed by protein family with protein family target IDs */
    private Map<String, SortedSet<Coupling>> couplingMap;

    // COMMAND-LINE OPTIONS

    /** output file from CouplesProcessor run containing couplings and scores */
    @Argument(index = 0, metaVar = "couples.tbl", usage = "output file containing coupling pairs and scores", required = true)
    private File inFile;

    /** directory of GTOs to update */
    @Argument(index = 1, metaVar = "genomeDir", usage = "directory of GTOs to update", required = true)
    private File gtoDir;

    @Override
    protected void setDefaults() {
        // Set the defaults for the base-class parameters.
        this.setDefaultConfiguration();
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Validate the base-class parameters.
        this.validateConfiguration();
        // Verify the input directory.
        if (! this.gtoDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.gtoDir + " is not found or invalid.");
        this.inputGenomes = new GenomeDirectory(this.gtoDir);
        log.info("Input directory is {}.", this.inputGenomes);
        // Verify the input file.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input couplings file " + this.inFile + " is not found or unreadable.");
        this.readCouplingMap();
        return true;
    }

    /**
     * Read the couplings from the input file.
     *
     * @throws IOException
     */
    private void readCouplingMap() throws IOException {
        try (TabbedLineReader inStream = new TabbedLineReader(this.inFile)) {
            // Find the key score columns.
            int sizeCol = inStream.findField("size");
            int strengthCol = inStream.findField("sim_distance");
            // Initialize the map.
            this.couplingMap = new HashMap<String, SortedSet<Coupling>>(10000);
            // Get the classifier for reading the family IDs.
            int count = 0;
            FeatureClass classifier = this.getClassifier();
            for (TabbedLineReader.Line line : inStream) {
                FeatureClass.Pair pair = classifier.readPair(line);
                int size = line.getInt(sizeCol);
                double strength = line.getDouble(strengthCol);
                this.addCoupling(pair.getClass1(), new Coupling(pair.getClass2(), size, strength));
                this.addCoupling(pair.getClass2(), new Coupling(pair.getClass1(), size, strength));
                count++;
            }
            log.info("{} couplings read from {}. {} families have couplings.", count, this.inFile,
                    this.couplingMap.size());
        }
    }

    /**
     * Add a coupling to the coupling map.
     *
     * @param classX	coupled class
     * @param coupling	coupling to that class
     */
    private void addCoupling(String classX, Coupling coupling) {
        SortedSet<Coupling> couplingList = this.couplingMap.computeIfAbsent(classX, k -> new TreeSet<Coupling>());
        couplingList.add(coupling);
    }

    @Override
    protected void runCommand() throws Exception {
        // Get the classifier and the neighbor-finder.
        FeatureClass classifier = this.getClassifier();
        NeighborFinder finder = this.getFinder();
        // Loop through the genomes.
        for (Genome genome : this.inputGenomes) {
            log.info("Processing genome {}.", genome);
            // Get the result list for this genome.
            List<FeatureClass.Result> gResults = classifier.getResults(genome);
            // Erase all the existing couplings.
            for (Feature feat : genome.getFeatures())
                feat.clearCouplings();
            // Loop through the results.
            int couplingCount = 0;
            int eligibleCount = 0;
            int isolatedCount = 0;
            int n = gResults.size() - 1;
            for (int i = 0; i < n; i++) {
                FeatureClass.Result resI = gResults.get(i);
                // Get all the classes coupled to this result's classes.
                Map<String, Coupling> couplingTable = new HashMap<String, Coupling>(10);
                for (String classI : resI) {
                    SortedSet<Coupling> coupledClasses = this.couplingMap.get(classI);
                    if (coupledClasses != null) {
                        for (Coupling coupling : coupledClasses)
                            couplingTable.put(coupling.getTarget(), coupling);
                    }
                }
                if (couplingTable.size() > 0) {
                    eligibleCount++;
                    // Now get all our neighbors.
                    Collection<FeatureClass.Result> neighbors = finder.getNeighbors(gResults, i);
                    if (neighbors.size() == 0) {
                        isolatedCount++;
                    } else {
                        Feature feat = genome.getFeature(resI.getFid());
                        for (FeatureClass.Result resJ : neighbors) {
                            // Find out if we have a coupling.  We keep the size and strength of the largest.
                            int size = 0;
                            double strength = 0.0;
                            for (String classJ : resJ) {
                                Coupling target = couplingTable.get(classJ);
                                if (target != null && target.getStrength() > strength) {
                                    size = target.getSize();
                                    strength = target.getStrength();
                                }
                            }
                            if (size > 0) {
                                feat.addCoupling(resJ.getFid(), size, strength);
                                Feature featJ = genome.getFeature(resJ.getFid());
                                featJ.addCoupling(resI.getFid(), size, strength);
                                couplingCount++;
                            }
                        }
                    }
                }
            }
            log.info("Genome {} had {} classifiable features, {} with eligible classes, with {} having no neighbors.",
                    genome, gResults.size(), eligibleCount, isolatedCount);
            if (couplingCount <= 0)
                log.warn("No couplings were found.");
            else {
                File outFile = this.inputGenomes.currFile();
                log.info("Updating to file {}. {} couplings were found.", outFile, couplingCount);
                genome.update(outFile);
            }
        }
        log.info("All done.");
    }

}
