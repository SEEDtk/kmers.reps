/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.IOException;
import java.io.PrintWriter;

import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command compares a list of universal roles (produced by UniRoleProcessor) to a directory of genomes.
 * It will output the protein families found for each role, plus a list of the roles missing from each of the
 * genomes.
 *
 * @author Bruce Parrello
 *
 */
public class MissingRoleProcessor extends BaseReportProcessor {

    @Override
    protected void setReporterDefaults() {
        // TODO code for setReporterDefaults

    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // TODO code for validateReporterParms

    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // TODO code for runReporter

    }

}
