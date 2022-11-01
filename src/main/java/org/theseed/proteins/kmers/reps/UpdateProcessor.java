/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.utils.ICommand;

/**
 * This script is used when an incremental update is being made to the genome database.  It takes as input
 * evaluation results for new genomes and folds them into existing results.  For a complete rerun of all
 * evaluation results, use the GenomeProcessor class instead.
 *
 * The positional parameters are the name of the output directory and the name of the original input directory.
 * The output directory must contain the updated "patric.eval.tbl" file containing the latest evaluation data.
 * The input directory must contain the repXX.ser files for the current representative genome databases.
 * The new genomes are processed to flesh out the repgen sets and then the updated files written to the
 * output directory.
 *
 * The following command-line options are supported.
 *
 * -v	show more detailed log messages
 * -b	batch size for PATRIC queries
 *
 * @author Bruce Parrello
 *
 */
public class UpdateProcessor extends BaseGenomeProcessor implements ICommand {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UpdateProcessor.class);
    /** input file for protein data factory */
    private File inFile;
    // COMMAND-LINE OPTIONS

    /** old evaluation directory */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory with existing genome files", required = true)
    private File inDir;

    @Override
    protected void setDefaults() {
        this.setBaseDefaults();
    }

    @Override
    protected boolean validateParms() throws IOException {
        this.checkParms();
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Verify that the output directory contains the input file.
        this.inFile = new File(this.getOutDir(), "patric.eval.tbl");
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Evaluation result file " + this.inFile + " is not found or unreadable.");
        // Verify that the input directory
        return true;
    }

    @Override
    public void runCommand() {
        try {
            // Read all the genomes from the input file and build protein data objects.
            initializeProteinData(inFile);
            // We need to create the FASTA files for the seed protein list and the
            // binning BLAST database.  We do that here.
            createFastaFiles();
            // Create a list of RepGenomeDb objects, one for each score limit.
            log.info("Reloading repgen sets.");
            ProteinDataFactory.restoreData(this, this.inDir);
            // Sort the genomes into repgen sets.
            collateGenomes();
            // Save all the repgen sets.
            saveRepGenSets();
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
            // Finally, save the repgen GTOs.
            // TODO check for GTO directoriesif (! this.noGTOs)
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void finishGenomeList(ProteinDataFactory gList) throws UnsupportedEncodingException, IOException {
        gList.finishList(this.getBatchSize());
    }

}
