package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.theseed.basic.BaseProcessor;
import org.theseed.genome.Feature;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceKmers;
import org.theseed.stats.QualityCountMap;
import org.slf4j.LoggerFactory;

/**
 * This is the primary class for processing a representative-genome database.  A FASTA file
 * is read from the standard input:  each record must have the key protein's feature ID as
 * the label, the genome name as the comment, and the key protein's sequence as the sequence.
 *
 * The positional parameter is the name of the representative-genome database's saved data file.
 * This is also a FASTA file, containing an empty-sequence header record and one record per
 * representative genome.
 *
 * A report will be produced to the standard output about the input genomes, specifying which are
 * represented, and the similarity score for the closest representative of each.  If the input
 * file is empty, no report will be produced.
 *
 * --create		create the representative-genome database from the specified file, which should have
 * 				the same format as the input file
 * --list		specifies a file to contain a listing of the ID and name of each representative genome
 * --null		the standard input will not be read
 * --verbose	progress messages will be produced on the standard error output
 * --stats		specifies a file to contain a report of the number of genomes found close to each representative
 *
 * When CREATE is specified, the following options are also used
 *
 * -m		minimum similarity threshold for representation (default 100)
 * -K		protein kmer size (default 8)
 * -p		key protein name (default "Phenylalanyl-tRNA synthetase alpha chain")
 * -x		maximum number of ambiguity characters allowed in a sequence; genomes with more than this number
 * 			are deferred to a second pass to reduce their likelihood of being selected (default 20)
 * --dna	if specified, the sequences are presumed to be DNA instead of amino acids
 *
 * @author Bruce Parrello
 *
 */
public class RepGenomeDbProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RepGenomeDbProcessor.class);
    /** number of genomes processed */
    private int genomesProcessed;
    /** start time of run */
    private long startTime;
    /** target representative-genome database */
    private RepGenomeDb repDB;
    /** saved time stamp */
    private long lastTime;

    // COMMAND LINE

    /** similarity threshold */
    @Option(name="-m", aliases={"--sim", "--minScore"}, metaVar="100", usage="similarity threshold for representation")
    private int threshold;

    /** kmer size */
    @Option(name="-K", aliases={"--kmer"}, metaVar="8", usage="protein kmer size")
    private void setKmer(int newSize) {
        ProteinKmers.setKmerSize(newSize);
    }

    /** key protein name */
    @Option(name="-p", aliases={"--prot", "--role"}, metaVar="\"role name\"", usage="name of the role for the key protein")
    private String protName;

    /** create command */
    @Option(name="--create", metaVar="fastaFile", usage="create the database from genomes in the specified file")
    private File createFile;

    /** list command */
    @Option(name="--list", metaVar="listFile", usage="list the representative genomes in the specified file")
    private File listFile;

    /** stats request */
    @Option(name="--stats", metaVar="statFile", usage="list statistics about genome representation on the specified file")
    private File statFile;

    /** input suppressed */
    @Option(name="--null", usage="do not read the standard input or produce a report")
    private boolean noInput;

    /** treat the sequences as DNA instead of amino acids */
    @Option(name="--dna", usage="input is DNA instead of AA")
    private boolean dnaMode;

    /** maximum number of ambiguity characters */
    @Option(name="-x", aliases={"--maxErr"}, metaVar="20", usage="maximum number of Ns (DNA) or Xs (amino acid)")
    private int maxErr;

    /** representative-genome file */
    @Argument(index=0, metaVar="repDbFile", usage="representative genome database file",
            required=true, multiValued=false)
    private File dbFile;

    /**
     * This method is called when we need to create the database.
     * @throws IOException
     */
    private void createRepDb() throws IOException {
        log.info("Creating rep db with K={} and threshold {} using {}.",
                ProteinKmers.kmerSize(), this.threshold, this.protName);
        this.repDB = new RepGenomeDb(this.threshold, this.protName);
        FastaInputStream inStream = new FastaInputStream(this.createFile);
        // We will put deferred sequences in here.
        ArrayList<Sequence> deferred = new ArrayList<Sequence>();
        // Compute the ambiguity character.
        char badChar = (this.dnaMode ? 'N' : 'X');
        // Loop through the input stream, processing representative genomes.  We can
        // do this in one statement, but we slow down to display progress.
        log.info("Processing input.");
        genomesProcessed = 0;
        this.lastTime = System.currentTimeMillis();
        startTime = lastTime;
        for (Sequence inSequence : inStream) {
            // Verify that this is a good sequence.
            if (StringUtils.countMatches(inSequence.getSequence(), badChar) >= this.maxErr) {
                deferred.add(inSequence);
            } else {
                processSequence(inSequence);
            }
        }
        if (deferred.size() > 0) {
           log.info("Processing {} deferred sequences.", deferred.size());
            for (Sequence inSequence : deferred) {
                processSequence(inSequence);
            }
        }
        if (log.isDebugEnabled()) showCreateProgress();
        inStream.close();
        // Save the database.
        log.info("Saving database to {}", dbFile);
        repDB.save(dbFile);
    }

    /**
     * Process a single sequence.
     *
     * @param inSequence	sequence to process
     */
    private void processSequence(Sequence inSequence) {
        RepGenome newGenome = new RepGenome(inSequence);
        repDB.checkGenome(newGenome);
        genomesProcessed++;
        if ((System.currentTimeMillis() - this.lastTime) > 10000) {
            // A minute since our last progress. Show more progress.
            if (log.isDebugEnabled()) showCreateProgress();
            this.lastTime = System.currentTimeMillis();
        }
    }

    /** Display our current progress during database creation. */
    private void showCreateProgress() {
        long rate = 0;
        if (genomesProcessed > 0) {
            rate = genomesProcessed * 1000 / (System.currentTimeMillis() - startTime);
        }
        log.debug("{} total genomes input, {} kept ({} genomes/second).",
                genomesProcessed, repDB.size(), rate);
    }

    @Override
    protected void setDefaults() {
        this.threshold = 100;
        this.protName = RepGenomeDb.DEFAULT_PROTEIN;
        this.dbFile = null;
        this.createFile = null;
        this.listFile = null;
        this.statFile = null;
        this.noInput = false;
        this.dnaMode = false;
        this.maxErr = 20;
    }

    @Override
    protected void validateParms() throws IOException {
        if (this.createFile == null && ! this.dbFile.exists())
            throw new FileNotFoundException("Create not specified and specified repDB file not found.");
    }

    @Override
    protected void runCommand() throws Exception {
        long initTime = System.currentTimeMillis();
        if (this.createFile != null) {
            // Here we must create the rep-genome database.
            createRepDb();
        } else {
            // Load the rep-genome database.
            log.info("Loading representative genome data from {}.", this.dbFile);
            this.repDB = RepGenomeDb.load(this.dbFile);
        }
        // Now we have our database in memory.  List the genomes, if necessary.
        if (this.listFile != null) {
            log.info("Writing representative genome names to {}.", this.listFile);
            PrintWriter listStream = new PrintWriter(listFile);
            SortedSet<RepGenome> allReps = repDB.sorted();
            listStream.println("genome_id\tname");
            for (RepGenome rep : allReps) {
                listStream.format("%s\t%s%n", rep.getGenomeId(), rep.getName());
            }
            listStream.close();
            log.info("{} genomes listed.", allReps.size());
        }
        // Only proceed if we have input data.
        if (! noInput) {
            // Now it's time to produce the report. Create an output header.
            log.info("Producing report.");
            System.out.println("genome_id\tname\trep_id\tscore");
            FastaInputStream inStream = new FastaInputStream(System.in);
            // Set up progress timing and stats.
            this.genomesProcessed = 0;
            int outliers = 0;
            QualityCountMap<RepGenome> repCounts = new QualityCountMap<RepGenome>();
            this.startTime = System.currentTimeMillis();
            // Loop through the input.
            for (Sequence inSeq : inStream) {
                // Test this genome.
                RepGenomeDb.Representation result = repDB.findClosest(inSeq);
                this.genomesProcessed++;
                // Classify the result.
                if (! result.isRepresented()) {
                    outliers++;
                    repCounts.setBad(result.getRepresentative());
                    log.debug("{} is an outlier with a score of {}.", inSeq.getLabel(), result.getSimilarity());
                } else {
                    repCounts.setGood(result.getRepresentative());
                }
                // Write it to the report.  Note we tweak the infinity value.
                int sim = result.getSimilarity();
                String simString = (sim == SequenceKmers.INFINITY ? "MAX" : Integer.toString(sim));
                String genomeId = Feature.genomeOf(inSeq.getLabel());
                System.out.format("%s\t%s\t%s\t%s%n", genomeId, inSeq.getComment(),
                        result.getRepresentative().getGenomeId(), simString);
                if (log.isInfoEnabled() && this.genomesProcessed % 100 == 0) {
                    long rate = 0;
                    if (genomesProcessed > 0) {
                        rate = this.genomesProcessed * 1000 /
                                (System.currentTimeMillis() - this.startTime);
                    }
                    log.info("{} genomes processed, {} outliers, {} genomes/second.",
                            this.genomesProcessed, outliers, rate);
                }
            }
            inStream.close();
            // End of report.
            if (this.statFile != null) {
                // Here the user wants statistics.
                log.info("Producing stat report on {}.", statFile.getPath());
                PrintWriter statStream = new PrintWriter(statFile);
                statStream.println("rep_id\trep_name\trepresented\toutliers");
                List<RepGenome> statReps = repCounts.bestKeys();
                for (RepGenome rep : statReps) {
                    statStream.format("%s\t%s\t%d\t%d%n", rep.getGenomeId(),
                            rep.getName(), repCounts.good(rep),
                            repCounts.bad(rep));
                }
                statStream.close();
                log.info("Statistics report done.");
            }
        }
        long totalDuration = (System.currentTimeMillis() - initTime) / 60000;
        log.info("Total duration is {} minutes.", totalDuration);
    }
}
