/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Coupling;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.coupling.FeatureClass.Pair;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.ParseFailureException;

/**
 * This command creates specially-augmented GTOs that can be used to drive the molecular machine processor.  Each
 * GTO will contain coupling information taken from the output of the CouplesProcessor with the SCORES option.
 *
 * The positional parameters are the name of the CouplesProcessor output file and the name of an input directory
 * containing GTOs.  The GTOs will be updated in place to contain coupling data.  In addition, a "couples.tbl"
 * file will be produced that lists all the instances of each coupling, and a "genomes.tbl" file will be produced
 * that lists all the genomes that contain couplings.  These files are used by the web applications.
 * For each coupling pair, the coupling file contains an output line consisting of the pair name and a
 * comma-delimited list of the participating feature pairs.  For each genome, the genome file contains an
 * output line consisting of the genome ID, the genome name, and the number of coupled features.
 *
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

    /**
     * This sorts coupling data sets from largest to smallest.
     */
    public class SetComparator implements Comparator<Entry<Pair, SortedSet<String>>> {

        @Override
        public int compare(Entry<Pair, SortedSet<String>> o1, Entry<Pair, SortedSet<String>> o2) {
            int retVal = o2.getValue().size() - o1.getValue().size();
            if (retVal == 0)
                retVal = o1.getKey().getClass1().compareTo(o2.getKey().getClass1());
            if (retVal == 0)
                retVal = o1.getKey().getClass2().compareTo(o2.getKey().getClass2());
            return retVal;
        }

    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PrepareProcessor.class);
    /** input directory */
    private GenomeDirectory inputGenomes;
    /** coupling map, keyed by protein family with protein family target IDs */
    private Map<String, SortedSet<Coupling>> couplingMap;
    /** map of coupled pairs to features */
    private Map<FeatureClass.Pair, SortedSet<String>> featureMap;


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
    protected boolean validateParms() throws IOException, ParseFailureException {
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
     * Read the couplings from the input file.  We fill in the coupling map and initialize all of
     * the feature map entries so that each coupling pair has an empty feature list.
     *
     * @throws IOException
     */
    private void readCouplingMap() throws IOException {
        try (TabbedLineReader inStream = new TabbedLineReader(this.inFile)) {
            // Find the key score columns.
            int sizeCol = inStream.findField("size");
            int strengthCol = inStream.findField("sim_distance");
            // Initialize the maps.
            this.couplingMap = new HashMap<String, SortedSet<Coupling>>(10000);
            this.featureMap = new HashMap<FeatureClass.Pair, SortedSet<String>>(10000);
            // Get the classifier for reading the family IDs.
            int count = 0;
            FeatureClass classifier = this.getClassifier();
            for (TabbedLineReader.Line line : inStream) {
                FeatureClass.Pair pair = classifier.readPair(line);
                int size = line.getInt(sizeCol);
                double strength = line.getDouble(strengthCol);
                this.addCoupling(pair.getClass1(), new Coupling(pair.getClass2(), size, strength));
                this.addCoupling(pair.getClass2(), new Coupling(pair.getClass1(), size, strength));
                this.featureMap.put(pair, new TreeSet<String>());
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
        // Open the genome output file.
        try (PrintWriter genomeStream = new PrintWriter(new File(this.gtoDir, "genomes.tbl"))) {
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
                                    if (target != null) {
                                        // Here we have a valid coupling.  Record it in the feature map.
                                        for (String classI : resI) {
                                            FeatureClass.Pair pair = classifier.new Pair(classI, classJ);
                                            SortedSet<String> coupledFeatures = this.featureMap.get(pair);
                                            String coupleString = resI.getFid() + ':' + resJ.getFid();
                                            coupledFeatures.add(coupleString);
                                        }
                                        // If it is the best coupling,
                                        if (target.getStrength() > strength) {
                                            size = target.getSize();
                                            strength = target.getStrength();
                                        }
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
                File outFile = this.inputGenomes.currFile();
                log.info("Updating to file {}. {} couplings were found.", outFile, couplingCount);
                genome.save(outFile);
                if (couplingCount <= 0)
                    log.warn("No couplings were found.");
                else {
                    genomeStream.format("%s\t%s\t%d%n", genome.getId(), genome.getName(), couplingCount);
                }
            }
            log.info("Sorting coupling-feature map.");
            List<Map.Entry<FeatureClass.Pair, SortedSet<String>>> featureMapList = new ArrayList<>(this.featureMap.entrySet());
            featureMapList.sort(new SetComparator());
            log.info("Unspooling coupling-feature map.");
            try (PrintWriter outStream = new PrintWriter(new File(this.gtoDir, "couples.tbl"))) {
                outStream.format("%s\tfeatures%n", classifier.getHeadings());
                for (Map.Entry<FeatureClass.Pair, SortedSet<String>> couplingData : featureMapList) {
                    SortedSet<String> fidSet = couplingData.getValue();
                    if (fidSet.size() > 0) {
                        String fidList = StringUtils.join(fidSet, ',');
                        outStream.format("%s\t%s%n", couplingData.getKey().toString(), fidList);
                    }
                }
            }
        }
        log.info("All done.");
    }

}
