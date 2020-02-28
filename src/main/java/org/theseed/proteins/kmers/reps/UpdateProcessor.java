/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.IOException;

import org.theseed.utils.ICommand;

/**
 * This script is used when an incremental update is being made to the genome database.  It takes as input
 * evaluation results for new genomes and folds them into existing results.  For a complete rerun of all
 * evaluation results, use the GenomeProcessor class instead.
 *
 * The positional parameters are the name of the output directory, the name of the input file, and the
 * name of the original input directory.  The input directory must contain the repXX.list.tbl files listing
 * the representatives of each genome, the PhenTrnaSyntAlph.fa file containing the seed protein DNA, and
 * the repXX.ser files containing the representative genome protein sequences.  This information is used
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

    @Override
    protected void setDefaults() {
        // TODO Auto-generated method stub

    }

    @Override
    protected boolean validateParms() throws IOException {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public void run() {
        // TODO Auto-generated method stub

    }

}
