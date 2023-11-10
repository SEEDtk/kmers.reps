/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.GenomeMultiDirectory;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

/**
 * This command creates a reference-genome database for evaluation (which is quite different from
 * one used for binning).  The database consists of a FASTA file of seed proteins and a master directory
 * of the genomes specified in the file.  The FASTA file is created by the subclasses of BaseGenomeProcessor.
 * For each genome whose seed protein is described in the FASTA file (the sequence ID will be the genome ID),
 * the genome is downloaded from PATRIC at the PROTEIN detail level.  This allows it to be fetched quickly
 * during evaluation.
 *
 * The positional parameters are the name of the input FASTA file and the output directory to contain
 * the reference-genome database.  The input FASTA file cannot already be in the output directory.
 * If the "--clear" option is not specified, then existing files already in the directory will be re-used
 * and obsolete files will be deleted.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --clear	erase the output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class BuildRefDbProcessor extends BaseProcessor {

    // FIELDS
    /** output master directory */
    private GenomeMultiDirectory outputGenomes;
    /** input genome list */
    private SortedMap<String, String> genomesIn;

    // COMMAND-LINE OPTIONS

    /** if specified, the output directory will be cleared before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** input FASTA file */
    @Argument(index = 0, metaVar = "refGenomes.fa", usage = "name of input FASTA file containing seed proteins with genome IDs for labels")
    private File inFile;

    /** output directory */
    @Argument(index = 1, metaVar = "outDir", usage = "output directory name")
    private File outDir;

    @Override
    protected void setDefaults() {
        this.clearFlag = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the input file.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " not found or unreadable.");
        // Get the absolute path of the output directory and verify it is not the input file's directory.
        String absOutDir = this.outDir.getAbsolutePath();
        if (this.inFile.getParentFile().getAbsolutePath().contentEquals(absOutDir))
            throw new IOException("Input file cannot be in output directory " + this.outDir + ".");
        // Now insure that the output directory exists and is ready.  This may fail if there is already
        // data in it.
        if (this.clearFlag || ! this.outDir.isDirectory()) {
            log.info("Initializing output directory {}.", this.outDir);
            this.outputGenomes = GenomeMultiDirectory.create(this.outDir, this.clearFlag);
        } else {
            log.info("Restarting in output directory {}.", this.outDir);
            this.outputGenomes = new GenomeMultiDirectory(this.outDir);
        }
        // Finally, we need to read in the genome list and copy the FASTA file to the master directory.
        log.info("Analyzing FASTA input file {}.", this.inFile);
        File fastaOutFile = new File(this.outDir, "refGenomes.fa");
        this.genomesIn = new TreeMap<String, String>();
        try (FastaInputStream fastaIn = new FastaInputStream(this.inFile);
                FastaOutputStream fastaOut = new FastaOutputStream(fastaOutFile)) {
            for (Sequence seq : fastaIn) {
                // Write the sequence to the output.
                fastaOut.write(seq);
                // Save the genome ID and name.
                this.genomesIn.put(seq.getLabel(), seq.getComment());
            }
        }
        log.info("{} reference genomes found in file.", this.genomesIn.size());
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // This part is very simple.  We need to get all the listed genomes from PATRIC and put them
        // in the output directory.  First, we need a PATRIC connection.
        log.info("Connecting to PATRIC.");
        var p3 = new P3Connection();
        // If the output directory has genomes in it, we need to do some pre-cleaning.  Note that we
        // do the check only for performance reasons.  The code would still work, it just wouldn't do
        // anything.
        if (this.outputGenomes.size() > 0) {
            // We must first remove obsolete genomes.
            List<String> oldIds = this.outputGenomes.getIDs().stream().filter(x -> ! this.genomesIn.containsKey(x))
                    .collect(Collectors.toList());
            log.info("{} obsolete genomes need to be removed.", oldIds.size());
            for (String oldId : oldIds)
                this.outputGenomes.remove(oldId);
            // Now we must reduce the genome list to contain only new genomes.
            var iter = this.genomesIn.entrySet().iterator();
            while (iter.hasNext()) {
                var genomeIn = iter.next();
                if (this.outputGenomes.contains(genomeIn.getKey()))
                    iter.remove();
            }
            log.info("{} genomes remain to be downloaded.", this.genomesIn.size());
        }
        // Now, we loop through the genome list.  The genomes were sorted by ID, so that the log
        // contains a sense of
        int processed = 0;
        int total = this.genomesIn.size();
        for (var genomeSpec : this.genomesIn.entrySet()) {
            String genomeId = genomeSpec.getKey();
            processed++;
            log.info("Downloading {} ({} of {}): {}.", genomeId, processed, total, genomeSpec.getValue());
            var genome = P3Genome.load(p3, genomeId, P3Genome.Details.PROTEINS);
            this.outputGenomes.add(genome);
        }
    }
}
