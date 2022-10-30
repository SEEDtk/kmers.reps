/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.theseed.utils.ICommand;

/**
 * This script is used when an incremental update is being made to the genome database.  It takes as input
 * a P3Eval directory as produced by the GenomeProcessor module and a new P3Eval directory containing
 * an updated "patric.eval.tbl".
 *
 * The positional parameters are the name of the output directory and the name of the original input
 * directory.  The output directory must contain the patric.eval.tbl file with the new evaluation results.
 * The input directory must contain the repXX.list.tbl files listing the representatives of each genome,
 * the PhenTrnaSyntAlph.fa file containing the seed protein DNA, the allProts.fa containing the seed
 * protein amino acid sequences, and the allSsu.fa containing the 16s SSU DNA sequences.  This information
 * is used to pre-load the ProteinData records and initialize the RepGenomeDb objects.
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
            // TODO perform the update
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
