/**
 *
 */
package org.theseed.genome.coupling;

import java.io.IOException;

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
public class DownloadProcessor extends BaseCouplingProcessor {

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
    protected void runCommand() throws Exception {
        // TODO Auto-generated method stub

    }

}
