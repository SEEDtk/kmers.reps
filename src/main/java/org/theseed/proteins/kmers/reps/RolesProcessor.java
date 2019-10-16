/**
 *
 */
package org.theseed.proteins.kmers.reps;

import org.theseed.utils.ICommand;

/**
 * This class processes a universal-role file against a representative-genome database.  The universal-role file is
 * tab-delimited, and for each genome it contains the genome ID, the seed protein sequence, and all of the genome's
 * singly-occurring roles.  For each representative genome, the output file contains the number of represented genomes
 * and a maximal set of such roles for which each genome contains at least a specified fraction of the roles in the set.
 *
 * The basic plan is to count the occurrences of each role and create an internal matrix of which roles each genome
 * contains.  We then create a set consisting of all the roles, and selectively remove the ones that occur least
 * frequently until the percentage limit is hit.
 *
 * @author Bruce Parrello
 *
 */
public class RolesProcessor implements ICommand {

    @Override
    public boolean parseCommand(String[] args) {
        // TODO parse roles processor
        return false;
    }

    @Override
    public void run() {
        // TODO execute roles processor

    }

}
