/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.utils.ICommand;
import org.theseed.utils.IntegerList;
import org.theseed.utils.ParseFailureException;

/**
 * This command processes the results from a quality run.  It expects as input the sorted results
 * file containing the lineage, the quality data, and the good-genome flag.  It will create
 * the master seed protein FASTA file, the master binning BLAST database, and the representative-genome
 * databases at level 10, 50, 100, and 200.
 *
 * The positional parameters are the name of the PATRIC master genome directory, the name of the output directory,
 * and the name of the input file.
 *
 * The following files are put into the output directory.
 *
 * The following command-line options are supported.
 *
 * -v	show more detailed log messages
 * -b	batch size for PATRIC queries
 *
 * --repSizes	comma-delimited list of score limits for each repGen set
 * --minLevel	minimum SSU rating level to keep (default SINGLE_SSU)
 * --noGtos		if specified, the GTO directories will be suppressed, greatly improving performance
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor extends BaseGenomeProcessor implements ICommand {

    /** list of RepGen set thresholds */
    private IntegerList repSizes;

    // COMMAND-LINE OPTIONS

    /** list of RepGen set scores */
    @Option(name = "--repSizes", aliases = { "--repScores", "--repSims" }, metaVar = "10,50,100,200",
            usage = "comma-delimited list of minimum scores for each RepGen set to create")
    private void setRepSizes(String repSizeString) {
        this.repSizes = new IntegerList(repSizeString);
    }

    /** if specified, no GTOs will be output */
    @Option(name = "--noGtos", aliases = { "--quick", "--noGTOs" }, usage = "if specified, no GTOs will be output")
    private boolean noGTOs;

    /** input tab-delimited file of genomes */
    @Argument(index = 2, metaVar = "inFile.tbl", usage = "input file (xxx.eval.tbl)", required = true)
    private File inFile;


    @Override
    protected void setDefaults() {
        this.setBaseDefaults();
        this.repSizes = new IntegerList(10, 50, 100, 200);
        this.noGTOs = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        this.checkParms();
        if (! this.inFile.canRead())
            throw new FileNotFoundException(this.inFile + " is not found or unreadable.");
        return true;
    }

    @Override
    public void runCommand() {
        try {
            // Read all the genomes from the input file and build protein data objects.
            initializeProteinData(this.inFile);
            // We need to create the FASTA files for the seed protein list and the
            // binning BLAST database.  We do that here.
            createFastaFiles();
            // Create a list of RepGenomeDb objects, one for each score limit.
            log.info("Creating repgen sets.");
            this.initRepGenSets(this.repSizes.size());
            for (int repSize : this.repSizes) {
                this.addRepGenSet(new RepGenomeDb(repSize, ProteinDataFactory.SEED_FUNCTION));
            }
            // Sort the genomes into repgen sets.
            collateGenomes();
            // Save all the repgen sets.
            saveRepGenSets();
            // Write out the protein Fasta file for the first set.  This is used to find
            // seed proteins.
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
            if (! this.noGTOs)
                saveGTOs();
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }


}
