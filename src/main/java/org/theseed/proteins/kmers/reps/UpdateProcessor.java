/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.kmers.ProteinDataFactory;

/**
 * This script is used when an incremental update is being made to the genome database.  It takes as input
 * evaluation results for new genomes and folds them into existing results.  For a complete rerun of all
 * evaluation results, use the GenomeProcessor class instead.
 *
 * The positional parameters are the name of the output directory and the name of the original input directory.
 * The output directory must contain the updated "patric.eval.tbl" file containing the latest evaluation data.
 * The input directory must contain the repXX.ser files for the current representative genome databases and
 * the "repXX.list.tbl" files containing the genome repgen placements from the previous evaluation run.
 * The new genomes are processed to flesh out the repgen sets and then the updated files written to the
 * output directory.
 *
 * The following command-line options are supported.
 *
 * -v	show more detailed log messages
 * -b	batch size for PATRIC queries
 * 
 * NOTE: due to a change in the way repgen sets are computed, this command no longer works
 *
 * @author Bruce Parrello
 *
 */
public class UpdateProcessor extends BaseGenomeProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(UpdateProcessor.class);
    /** input file for protein data factory */
    private File inFile;
    // COMMAND-LINE OPTIONS

    /** old evaluation results directory */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory with existing genome files", required = true)
    private File inDir;

    @Override
    protected void setDefaults() {
        this.setBaseDefaults();
    }

    @Override
    protected void validateParms() throws IOException {
        this.checkParms();
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Verify that the output directory contains the input file.
        this.inFile = new File(this.getOutDir(), "patric.eval.tbl");
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Evaluation result file " + this.inFile + " is not found or unreadable.");
        // Verify that the input directory
    }

    @Override
    public void runCommand() {
        try {
            if (! this.isFullRegenRequired()) {
                // Here we need to reload the repgen sets from the output directory so we can download
                // the GTOs.
                log.info("Reloading new repgen sets.");
                ProteinDataFactory.restoreData(this, this.getOutDir());
            } else {
                // Read all the genomes from the input file and build protein data objects.
                initializeProteinData(inFile);
                // We need to create the FASTA files for the seed protein list and the
                // binning BLAST database.  We do that here.
                createFastaFiles();
                // Create a list of RepGenomeDb objects, one for each score limit.
                log.info("Reloading old repgen sets.");
                ProteinDataFactory.restoreData(this, this.inDir);
                // Remove any representative genome that doesn't exist.
                int deleteCount = 0;
                var repGenSets = this.getRepGenSets();
                for (RepGenomeDb repDb : repGenSets) {
                    // Get an iterator so we can remove the missing ones as we scan.
                    var iter = repDb.iterator();
                    while (iter.hasNext()) {
                        RepGenome rep = iter.next();
                        if (! this.isPresent(rep.getGenomeId())) {
                            deleteCount++;
                            log.info("Deleting {} from {}.", rep.toString(), repDb.toString());
                            iter.remove();
                        }
                    }
                }
                log.info("{} representatives had to be deleted.", deleteCount);
                // Save the repgen set levels in the protein data factory.
                saveRepLevels();
                // Determine which genomes have valid placements in all the repgen sets.
                Set<String> skipSet = this.computeSkipSet(repGenSets);
                // Sort the genomes into repgen sets.
                collateGenomes(skipSet);
                // Save all the repgen sets.
                saveRepGenSets();
                // Write the protein data report.
                writeGenomeReport();
                // Write out the protein Fasta file for the first set.  This is used to find
                // seed proteins during binning.
                writeSeedProt();
                // Assign genomes to repgen sets.
                log.info("Assigning genomes to repgen sets.");
                writeListFiles();
                // Now we write the protein FASTA files and the stats files.
                writeProteinFasta();
                // Next, the reference-genome finder for evaluation.
                writeRefGenomeFasta();
                // Now we produce the repFinder file used to find close genomes.
                writeRepFinder();
                // Write the good-genome list.
                writeGoodGenomes();
            }
            // Finally, save the repgen GTOs.
            if (this.isGTOsRequested()) {
                log.info("Creating repgen GTO directories.");
                setupGTOs();
                log.info("Copying old GTOs and identifying new ones.");
                copyGTOs(this.inDir);
                log.info("Downloading new GTOs.");
                saveNewGTOs();
            }
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * This method runs through the repgen set "list.repXX.tbl" files in the input directory,
     * to determine which are already placed.  A genome is placed if its representative has not
     * been deleted.
     *
     * @param repGenSets	list of repgen sets of interest
     *
     * @return a set of the IDs for the genomes that do not need to be checked when computing new representatives
     *
     * @throws IOException
     */
    private Set<String> computeSkipSet(List<RepGenomeDb> repGenSets) throws IOException {
        // This will count the number of repgen sets where each genome has been placed.
        CountMap<String> counts = new CountMap<String>();
        for (RepGenomeDb repdb : repGenSets) {
            // Get the list file for this repgen set.
            File inFile = new File(this.inDir, repdb.getListFileName());
            if (! inFile.exists())
                log.warn("WARNING: List file for {} not found in {}.", repdb, this.inDir);
            else {
                // The file exists, so we can analyze it.
                log.info("Analyzing list file {}.", inFile);
                try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
                    int gCol = inStream.findField("genome_id");
                    int rCol = inStream.findField("rep_id");
                    for (var line : inStream) {
                        String genomeId = line.get(gCol);
                        if (this.isPresent(genomeId)) {
                            String repId = line.get(rCol);
                            if (this.isPresent(repId))
                                counts.count(genomeId);
                        }
                    }
                }
            }
        }
        // Now we return every genome with a correct count.
        final int setCount = repGenSets.size();
        Set<String> retVal = counts.counts().stream().filter(x -> x.getCount() >= setCount)
                .map(x -> x.getKey()).collect(Collectors.toSet());
        log.info("{} genomes eligible to be skipped.", retVal.size());
        return retVal;
    }

    @Override
    protected void finishGenomeList(ProteinDataFactory gList) throws UnsupportedEncodingException, IOException {
        gList.finishList(this.getBatchSize());
    }

}
