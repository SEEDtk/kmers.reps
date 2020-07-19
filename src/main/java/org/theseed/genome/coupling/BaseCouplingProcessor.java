package org.theseed.genome.coupling;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.utils.BaseProcessor;

/**
 * This is a base class for coupling-related commands.  It provides common methods for manipulating feature classes
 * and neighbors.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseCouplingProcessor extends BaseProcessor {

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

    /**
     * Set the defaults for the common parameters.
     */
    protected void setDefaultConfiguration() {
        this.classType = FeatureClass.Type.PGFAMS;
        this.neighborType = NeighborFinder.Type.ADJACENT;
        this.maxGap = 5000;
    }

    /**
     * Validate the common parameters.
     */
    protected void validateConfiguration() {
        if (this.maxGap < 1)
            throw new IllegalArgumentException("Invalid maximum gap size.  Must be at least 1.");
        // Create the feature classifier.
        this.classifier = this.classType.create();
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

}
