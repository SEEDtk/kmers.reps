/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;

/**
 * This class reads a representative-genome database and computes the distance between each pair of genomes
 * in the database.  The positional parameter is the name of the representative-genome database file.  The
 * output will be a three-column table on the standard output.
 *
 * The command-line parameters are as follows.
 *
 * -h 	display command-line usage
 * -v	display progress messages on log
 *
 * @author Bruce Parrello
 *
 */
public class DistanceMatrixProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DistanceMatrixProcessor.class);
    /** target representative-genome database */
    private RepGenomeDb repDB;

    // COMMAND-LINE OPTIONS

    /** representative-genome file */
    @Argument(index=0, metaVar="repDbFile", usage="representative genome database file",
            required=true, multiValued=false)
    private File dbFile;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.dbFile.canRead())
            throw new FileNotFoundException("RepDB file " + this.dbFile + " not found or unreadable.");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Load the rep-genome database.
        log.info("Loading representative genome data from {}.", this.dbFile);
        this.repDB = RepGenomeDb.load(this.dbFile);
        // Write the output header.
        System.out.println("genome_1\tgenome_2\tdistance\tsim");
        // Compute the number of pairs to process.
        int totalPairs = this.repDB.size() * (this.repDB.size() + 1) / 2;
        log.info("{} total genomes, requiring {} pairs to be processed.", this.repDB.size(), totalPairs);
        // Loop through the representatives.
        long startTime = System.currentTimeMillis();
        int processed = 0;
        double distSum = 0;
        double distSqSum = 0;
        double simSum = 0;
        double simSqSum = 0;
        int zeroes = 0;
        for (RepGenome rep : this.repDB) {
            for (RepGenome rep2 : this.repDB) {
                // Compare the two genomes to insure each pair is only processed one.
                if (rep.compareTo(rep2) < 0) {
                    double distance = rep.distance(rep2);
                    int sim = rep.similarity(rep2);
                    distSum += distance;
                    distSqSum += distance * distance;
                    simSum += sim;
                    simSqSum += sim * sim;
                    if (sim == 0) zeroes++;
                    System.out.format("%s\t%s\t%8.6f\t%d%n", rep.getGenomeId(), rep2.getGenomeId(), distance, sim);
                    processed++;
                    if (log.isInfoEnabled() && processed % 50000 == 0) {
                        double rate = (double) processed / ((System.currentTimeMillis() - startTime) / 1000.0);
                        log.info("{} of {} pairs processed, {} pairs/second.", processed, totalPairs,
                                rate);
                    }
                }
            }
        }
        long duration = (System.currentTimeMillis() - startTime) / 1000;
        log.info("{} pairs processed in {} seconds.", processed, duration);
        double distMean = distSum / processed;
        double distSTdev = Math.sqrt(distSqSum / processed - distMean * distMean);
        double simMean = simSum / processed;
        double simSTdev = Math.sqrt(simSqSum / processed - simMean * simMean);
        log.info("Mean distance is {} with standard deviation {}.", distMean, distSTdev);
        log.info("Mean similarity is {} with standard deviation {}.", simMean, simSTdev);
        log.info("{} pairs had 0 similarity.", zeroes);
    }

}
