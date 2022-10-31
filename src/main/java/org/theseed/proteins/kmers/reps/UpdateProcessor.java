/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

import org.kohsuke.args4j.Argument;
import org.theseed.utils.ICommand;

/**
 * This script is used when an incremental update is being made to the genome database.  It takes as input
 * evaluation results for new genomes and folds them into existing results.  For a complete rerun of all
 * evaluation results, use the GenomeProcessor class instead.
 *
 * The positional parameters are the name of the output directory and the name of the original input directory.
 * The output directory must contain the updated "patric.eval.tbl" file containing the latest evaluation data.
 * The input directory must contain the repXX.list.tbl files listing the representatives of each genome,
 * the PhenTrnaSyntAlph.fa file containing the seed protein DNA, the allProts.fa file containing the seed
 * protein sequences, and the allSsu.fa file containing the SSU DNA sequences.  This information is used
 * to preload the ProteinData records and initialize the RepGenomeDb objects.  The new genomes are then
 * read in and the updated files written to the output directory.
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
        return true;
    }

    @Override
    public void runCommand() {
        try {
            // TODO update old P3Eval into new one
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void finishGenomeList(ProteinDataFactory gList) throws UnsupportedEncodingException, IOException {
        // TODO code for finishGenomeList

    }

}
