/**
 *
 */
package org.theseed.proteins.kmers.reps;

import org.theseed.utils.ICommand;

/**
 * This class processes a universal-role file to create a rep-genome based completeness engine.
 * For each group, we need the group ID, the similarity score (0 for taxonomy-based groups), the
 * seed protein (empty for taxonomy-based groups), and a comma-delimited list of the universal roles.
 * Together with a file mapping role names to IDs, this is enough to build a completeness processor.
 *
 * The input file consists of one or more groups.  Each group has a header, one record per genome,
 * and a trailer.  Every record is tab-delimited.  The header contains the group ID, score, and
 * seed protein sequence.  Each data record contains a genome ID and a comma-delimited role list.
 * The trailer is simply the string "//".  Note that none of the fields contain spaces, so a simple
 * java Scanner class will parse the input correctly.
 *
 * The command-line options are as follows.
 *
 * -h	display the parameters
 * -v	show progress on STDERR
 * -i	specifies a file name to use for the standard input
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
