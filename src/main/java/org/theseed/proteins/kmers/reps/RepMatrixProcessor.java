/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Option;
import org.theseed.basic.BaseProcessor;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * This command reads a FASTA file and sorts them by distance.  It then produces a list
 * of sequences such that every sequence is no more than a specified maximum distance
 * from at least one sequence in the list.  The result is a set of sequences that can
 * be used to determine if other sequences belong in the list.
 *
 * The command-line parameters are as follows.
 *
 * -K		protein kmer size (default 8)
 * -v		display progress messages on STDERR
 * -i		input FASTA file (default is STDIN)
 * -m		maximum acceptable distance
 *
 * @author Bruce Parrello
 *
 */
public class RepMatrixProcessor extends BaseProcessor {

    // FIELDS

    /** input FASTA stream */
    FastaInputStream inStream;
    /** output FASTA stream */
    FastaOutputStream outStream;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name="-K", aliases={"--kmer"}, metaVar="8", usage="protein kmer size (default 8)")
    private void setKmer(int newSize) {
        ProteinKmers.setKmerSize(newSize);
    }

    /** input file (if not using STDIN) */
    @Option(name = "--input", aliases = { "-i" }, usage = "input file (if not STDIN)")
    private File inFile;

    /** maximum acceptable distance */
    @Option(name = "-m", aliases = { "--max", "--maxDist" }, usage = "maximum distance for a group")
    private double maxDist;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.maxDist = 0.995;
        this.setKmer(8);
    }

    @Override
    protected void validateParms() throws IOException {
        // Set up the input FASTA.
        if (this.inFile != null)
            this.inStream = new FastaInputStream(this.inFile);
        else
            this.inStream = new FastaInputStream(System.in);
        // Set up to write the FASTA.
        this.outStream = new FastaOutputStream(System.out);
    }

    @Override
    protected void runCommand() throws Exception {
        try {
            // The sequences will be stored in here.
            List<SequenceInfo> sequenceList = new ArrayList<SequenceInfo>(1000);
            // Loop through the input.
            for (Sequence seq : inStream) {
                log.debug("Processing sequence {}.", seq.getLabel());
                // Create a sequence-info object for this sequence and compare it to all previous sequences.
                SequenceInfo info = new SequenceInfo(seq);
                for (SequenceInfo other : sequenceList)
                    info.storeComparison(other);
                sequenceList.add(info);
            }
            // Sort by the mean distance.
            log.info("Sorting distances.");
            List<SequenceInfo> sorted  = sequenceList.stream().sorted().collect(Collectors.toList());
            // Divide the sequences into groups.  We will store the saved sequences in this
            // list.
            sequenceList.clear();
            for (SequenceInfo info : sorted) {
                // If this is a singleton, write it out.
                if (info.getMin() > this.maxDist) {
                    log.debug("{} is a singleton group.", info.getId());
                } else {
                    // Find out if we are close to any existing sequence.
                    SequenceInfo found = null;
                    Iterator<SequenceInfo> iter = sequenceList.iterator();
                    while (iter.hasNext() && found == null) {
                        SequenceInfo other = iter.next();
                        if (info.getDistance(other) <= this.maxDist) found = other;
                    }
                    if (found == null) {
                        // This sequence starts a new group.
                        sequenceList.add(info);
                        outStream.write(info.getSeq());
                        log.debug("{} starts a new group.", info.getId());
                    } else {
                        // This sequence is already in a group.
                        log.debug("{} is close to {}.", info.getId(), found.getId());
                    }
                }
            }
            System.err.println(sequenceList.size() + " sequences form the representation set.");
        } finally {
            inStream.close();
            outStream.close();
        }
    }

}
