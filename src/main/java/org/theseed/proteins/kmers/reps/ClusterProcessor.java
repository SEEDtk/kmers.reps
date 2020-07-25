/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.IOException;

import org.theseed.utils.BaseProcessor;

/**
 * This command runs through the output of CouplesProcessor and tries to organize the couplings into groups.  It takes as input the SCORES
 * output from the original couplings and remembers the size and coupling probabilities for each.  It will attempt to form groups where the
 * chance of seeing all the functions together is higher than a certain threshold.
 *
 * For any particular group, we track the total probability of the group.  When we add a new member to the group, we multiply its probability
 * of occurring with the
 *
 * @author Bruce Parrello
 *
 */
public class ClusterProcessor extends BaseProcessor {

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
