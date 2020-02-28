/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import org.kohsuke.args4j.Option;
import org.theseed.proteins.ProteinDataFactory;
import org.theseed.utils.ICommand;
import org.theseed.utils.IntegerList;

/**
 * This command processes the results from a quality run.  It expects as input the sorted results
 * file containing the lineage, the quality data, and the good-genome flag.  It will create
 * the master seed protein FASTA file, the master binning BLAST database, and the representative-genome
 * databases at level 10, 50, 100, and 200.
 *
 * The positional parameters are the name of the output directory and the name of the input file.
 *
 * The following files are put into the output directory.
 *
 * The following command-line options are supported.
 *
 * -v	show more detailed log messages
 * -b	batch size for PATRIC queries
 *
 * --repSizes	comma-delimited list of score limits for each repGen set
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor extends BaseGenomeProcessor implements ICommand {

    /** list of RepGen set thresholds */
    private IntegerList repSizes;
    

    // COMMAND-LINE OPTIONS

    /** list of RepGen set scores */
    @Option(name = "--repSizes", aliases = { "--repScores", "--repSims" }, metaVar = "50,100,200,250",
            usage = "comma-delimited list of minimum scores for each RepGen set to create")
    private void setRepSizes(String repSizeString) {
        this.repSizes = new IntegerList(repSizeString);
    }

    @Override
    protected void setDefaults() {
        this.batchSize = 500;
        this.repSizes = new IntegerList(10, 50, 100, 200);
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.outDir.exists()) {
            // Insure we have an output directory.
            log.info("Creating output directory {}.", this.outDir);
            if (! this.outDir.mkdir())
                throw new IOException("Could not create output directory " + this.outDir);
        } else if (! this.outDir.isDirectory()) {
            throw new FileNotFoundException("Invalid output directory " + this.outDir);
        }
        if (! this.inFile.canRead())
            throw new FileNotFoundException(this.inFile + " is not found or unreadable.");
        return true;
    }

    @Override
    public void run() {
        try {
            // Read all the genomes from the input file.
            initializeProteinData();
            // We need to create the FASTA files for the seed protein list and the
            // binning BLAST database.  We do that here.
            createFastaFiles();
            // Create a list of RepGenomeDb objects, one for each score limit.
            log.info("Creating repgen sets.");
            this.repGenSets = new ArrayList<RepGenomeDb>(this.repSizes.size());
            for (int repSize : this.repSizes) {
                this.repGenSets.add(new RepGenomeDb(repSize, ProteinDataFactory.SEED_FUNCTION));
            }
            // Sort the genomes into repgen sets.
            sortGenomes();
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
            // Now we produce the repFinder file used to find close genomes.
            writeRepFinder();
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
