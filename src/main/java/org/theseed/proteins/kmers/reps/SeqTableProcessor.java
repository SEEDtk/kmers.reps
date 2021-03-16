/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.locations.Location;
import org.theseed.proteins.Function;
import org.theseed.utils.ParseFailureException;
import org.theseed.utils.RestartableBaseProcessor;

/**
 * This command reads a directory of GTOs and outputs a table of genome ID, seed protein sequence, and SSU rRNA
 * nucleotide sequence.  These can be used by other programs to find the closest representative of an input
 * genome fragment.
 *
 * The positional parameters are the name of the input genome directories.  They must contain full GTO files.  The
 * output table will be written to the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed log messages
 *
 * --resume		if specified, the name of an output file; it is presumed the process was interrupted and it will be
 * 				resumed
 *
 * @author Bruce Parrello
 *
 */
public class SeqTableProcessor extends RestartableBaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SeqTableProcessor.class);
    /** match pattern for SSU rRNA */
    public static final Pattern SSU_R_RNA = Pattern.compile("SSU\\s+rRNA|Small\\s+Subunit\\s+(?:Ribosomal\\s+r)?RNA|ssuRNA|16S\\s+(?:r(?:ibosomal\\s+)?)?RNA", Pattern.CASE_INSENSITIVE);
    /** match pattern for seed protein */
    public static final Pattern SEED_PROTEIN = Pattern.compile("Phenylalanyl-tRNA\\s+synthetase\\s+alpha\\s+chain(?:\\s+\\(E[^)]+\\))?", Pattern.CASE_INSENSITIVE);
    /** list of sequence names */
    private static final String[] funSeqNames = new String[] { "seed_protein", "ssu_rna" };

    // COMMAND-LINE OPTIONS

    /** list of genome directories to process */
    @Argument(index = 0, metaVar = "gtoDir1 gtoDir2 ...", usage = "GTO directories to scan")
    private List<File> gtoDirs;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify each of the genome directories.
        for (File gtoDir : gtoDirs) {
            if (! gtoDir.isDirectory())
                throw new FileNotFoundException("Input directory " + gtoDir + " is not found or invalid.");
        }
        // Setup the output/resume file.
        String header = "genome_id\t" + StringUtils.join(funSeqNames, "\t");
        this.setup(header, "genome_id");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        int badCount = 0;
        int totCount = 0;
        for (File gtoDir : this.gtoDirs) {
            log.info("Processing directory {}.", gtoDir);
            GenomeDirectory genomes = new GenomeDirectory(gtoDir);
            for (Genome genome : genomes) {
                String genomeId = genome.getId();
                if (! this.isProcessed(genomeId)) {
                    String[] seqsFound = this.processGenome(genome);
                    // Verify that we found everything.  We do this the slow way to provide more useful
                    // error messages.
                    boolean ok = true;
                    for (int i = 0; i < funSeqNames.length; i++) {
                        if (seqsFound[i].isEmpty()) {
                            ok = false;
                            log.warn("WARNING: {} is missing a {} sequence.", genome, funSeqNames[i]);
                        }
                    }
                    if (ok)
                        this.println(genomeId + "\t" + StringUtils.join(seqsFound, "\t"));
                    else
                        badCount++;
                    this.markProcessed(genomeId);
                }
                totCount++;
            }
        }
        log.info("{} of {} genomes were missing sequences.", badCount, totCount);
    }

    /**
     * Search this genome for a seed protein and an SSU RNA, and write them to the output.  We keep the
     * longest sequence in each case.
     *
     * @param genome	genome to search.
     *
     * @return an array of the sequences found
     */
    private String[] processGenome(Genome genome) {
        // Initalize the output.
        String[] retVal = new String[funSeqNames.length];
        Arrays.fill(retVal, "");
        // Loop through the features.
        log.info("Searching {}.", genome);
        for (Feature feat : genome.getFeatures()) {
            String function = Function.commentFree(feat.getPegFunction());
            switch (feat.getType()) {
            case "CDS" :
                // For a peg, we check for the seed protein.
                if (SEED_PROTEIN.matcher(function).matches()) {
                    String proposed = feat.getPegFunction();
                    if (proposed.length() > retVal[0].length())
                        retVal[0] = proposed;
                }
                break;
            case "rna" :
                // For an RNA, we check for an SSU.
                if (SSU_R_RNA.matcher(function).find()) {
                    Location loc = feat.getLocation();
                    if (loc.getLength() > retVal[1].length())
                        retVal[1] = genome.getDna(feat.getLocation());
                }
                break;
            }
        }
        return retVal;
    }

}
