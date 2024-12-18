package org.theseed.genome.coupling;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;

/**
 * This is a base class for coupling-related commands.  It provides common methods for manipulating feature classes
 * and neighbors.
 *
 * The following command-line options are supported.
 *
 *  -h	display command-line usage
 *  -v	display more frequent log messages
 *  -t	type of feature classification used in the input (default PGFAMS)
 *  -n	algorithm for determining the feature neighborhood (default CLOSE)
 *  -d	maximum acceptable distance for features to be considered neighbors (default 5000)
 *
 *  --famFile	name of a file containing protein family definitions (for FILE_FAMILY or LOCAL_TEST classifier)
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseCouplingProcessor extends BaseProcessor implements FeatureClass.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BaseCouplingProcessor.class);
    /** feature classifier */
    private FeatureClass classifier;
    /** neighbor finder */
    private NeighborFinder finder;

    // COMMAND-LINE OPTIONS

    /** maximum gap distance */
    @Option(name = "-d", aliases = { "--distance", "--gap" }, metaVar = "4500",
            usage = "maximum acceptable gap between neighboring features")
    private int maxGap;

    /** classification method */
    @Option(name = "-t", aliases = { "--type" }, usage = "type of classification for features")
    private FeatureClass.Type classType;

    /** algorithm for finding neighboring features */
    @Option(name = "-n", aliases = { "--neighbors" }, usage = "algorithm for finding neighboring features")
    private NeighborFinder.Type neighborType;

    /** protein-family definition file */
    @Option(name = "--famFile", usage = "optional protein-family definition file")
    private File famDefinitionFile;

    /**
     * Set the defaults for the common parameters.
     */
    protected void setDefaultConfiguration() {
        this.classType = FeatureClass.Type.PGFAMS;
        this.neighborType = NeighborFinder.Type.CLOSE;
        this.maxGap = 5000;
        this.famDefinitionFile = null;
    }

    /**
     * Validate the common parameters.
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    protected void validateConfiguration() throws ParseFailureException, IOException {
        if (this.maxGap < 1)
            throw new ParseFailureException("Invalid maximum gap size.  Must be at least 1.");
        // Verify the family definition file.
        if (this.famDefinitionFile != null && ! this.famDefinitionFile.canRead())
            throw new FileNotFoundException("Family definintion file " + this.famDefinitionFile + " is not found or unreadable.");
        // Create the feature classifier.
        this.classifier = this.classType.create(this);
        // Create the neighbor finder.
        this.finder = this.neighborType.create(this);
        log.info("Feature classification is {}. Neighborhood algorithm is {}.", this.classType, this.neighborType);
    }

    /**
     * @return the feature-classification object for this processor
     */
    public FeatureClass getClassifier() {
        return this.classifier;
    }

    /**
     * @return the neighbor-finder object for this processor
     */
    public NeighborFinder getFinder() {
        return this.finder;
    }

    /**
     * @return the maximum gap distance
     */
    public int getMaxGap() {
        return this.maxGap;
    }

    /**
     * @return the name of the file containing the protein family definitions
     */
    public File getFamFile() {
        return this.famDefinitionFile;
    }

}
