/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * This command will read through the finder file for a SOUR, and output genomes that are
 * at various distributed distances from the first genome. The total distance range (0 to 1)
 * will be divided into N different buckets, and at most M genomes will be put in each bucket,
 * where N and M are command-line options. The full genome list will be output on the standard
 * output. Note that we will deliberately ignore distances of 1.0, which is most of them.
 *
 * This method will generally select a small subset of the total list of genomes. To vary the
 * subset chosen, you can specify a different starting genome. The default starting genome is
 * the first one in the finder file. To specify a different default, indicate a file containing
 * just the starting genome's entry from the finder file.
 *
 * The positional parameter is the name of the finder file for the SOUR in question.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for genome list (if not STDOUT)
 * -K	kmer size for distance (default 20)
 *
 * --buckets	number of buckets in which to divide the distance range (default 50)
 * --size		maximum number of genomes to allow in a bucket (default 4)
 * --primer		optional file containing the starting genome finder entry
 *
 * @author Bruce Parrello
 *
 */
public class DistributedDistanceProcessor extends BaseReportProcessor {

	// FIELDS
	/** logging facility */
	protected static Logger log = LoggerFactory.getLogger(DistributedDistanceProcessor.class);
	/** kmers for start protein */
	private DnaKmers sourceKmers;
	/** array of upper bucket limits */
	private double[] bucketMax;
	/** feature ID of the primer protein */
	private String primerID;
	/** list of genome buckets; each bucket is a map from genome ID to name */
	private List<Map<String, String>> bucketList;
	/** number of full buckets */
	private int fullCount;
	/** number of genomes kept */
	private int keptCount;
	/** number of empty buckets */
	private int emptyCount;
	/** number of unit distances */
	private int farCount;

	// COMMAND-LINE OPTIONS

	/** DNA kmers size */
	@Option(name = "--kmerSize", aliases = { "-K" }, metaVar = "23", usage = "kmer size for distance computation")
	private int kmerSize;

	/** number of distance buckets to use for genomes */
	@Option(name = "--buckets", metaVar = "100", usage = "number of distance category buckets to use")
	private int numBuckets;

	/** maximum bucket size */
	@Option(name = "--size", metaVar = "12", usage = "maximum bucket size")
	private int bucketSize;

	/** optional file containing finder entry for primer DNA sequence */
	@Option(name = "--primer", metaVar = "primer.fna", usage = "if specified, a file containing the primer protein finder entry")
	private File primerFasta;

	/** finder file name */
	@Argument(index = 0, metaVar = "finderProtFile.fna", usage = "finder file for SOUR protein", required = true)
	private File finderFasta;

	@Override
	protected void setReporterDefaults() {
		this.kmerSize = 20;
		this.numBuckets = 50;
		this.bucketSize = 4;
		this.primerFasta = null;
	}

	@Override
	protected void validateReporterParms() throws IOException, ParseFailureException {
		// Insure the finder file exists.
		if (! this.finderFasta.canRead())
			throw new FileNotFoundException("Finder file " + this.finderFasta + " is not found or unreadable.");
		// If no primer is specified, use the finder file itself.
		if (this.primerFasta == null)
			this.primerFasta = finderFasta;
		else if (! this.primerFasta.canRead())
			throw new FileNotFoundException("Primer file " + this.primerFasta + " is not found or unreadable.");
		// Validate the numeric parameters.
		if (this.kmerSize < 2)
			throw new ParseFailureException("Kmer size must be at least 2.");
		if (this.numBuckets < 1)
			throw new ParseFailureException("Number of buckets must be at least 1.");
		if (this.bucketSize < 1)
			throw new ParseFailureException("Bucket size must be at least 1.");
		// Set the kmer size.
		DnaKmers.setKmerSize(this.kmerSize);
		log.info("Kmer size is {}.", this.kmerSize);
		// Create the output bucket list.
		int hashSize = this.bucketSize * 5 / 3 + 1;
		this.bucketList = IntStream.range(0, this.numBuckets).mapToObj(i -> new HashMap<String, String>(hashSize))
				.collect(Collectors.toList());
		// Create the array of bucket maximums. Note we fix the last entry to prevent roundoff error.
		double width = 1.0 / this.numBuckets;
		this.bucketMax = IntStream.range(1, this.numBuckets+1).mapToDouble(i -> i*width).toArray();
		this.bucketMax[this.numBuckets - 1] = 1.0;
		// Denote every bucket is empty.
		this.fullCount = 0;
		this.emptyCount = this.numBuckets;
		this.keptCount = 0;
		this.farCount = 0;
	}

	@Override
	protected void runReporter(PrintWriter writer) throws Exception {
		int inCount = 0;
		// Open the primer FASTA file and get the first sequence.
		log.info("Reading primer protein from {}.", this.primerFasta);
		try (FastaInputStream inStream = new FastaInputStream(this.primerFasta)) {
			if (! inStream.hasNext())
				throw new IOException("Primer FASTA file " + this.primerFasta + " is empty.");
			this.readFirstSequence(inStream);
			inCount++;
		}
		// Open the input FASTA file.
		try (FastaInputStream inStream = new FastaInputStream(this.finderFasta)) {
			// Set up some progress indicator stuff for logging. We will output progress every
			// 10 seconds.
			long lastMsg = System.currentTimeMillis();
			// Now loop through the other sequences, adding them to buckets until all the buckets
			// fill or we run out of sequences.
			while (inStream.hasNext() && this.fullCount < this.numBuckets) {
				Sequence seq = inStream.next();
				// Only proceed if this is not the primer.
				if (! seq.getLabel().equals(this.primerID)) {
					DnaKmers targetKmers = new DnaKmers(seq.getSequence());
					double distance = this.sourceKmers.distance(targetKmers);
					// Find the bucket for this distance.
					OptionalInt iBucket = IntStream.range(0, this.numBuckets)
							.filter(i -> this.bucketMax[i] > distance).findFirst();
					// If there is no bucket, record a far genome. Otherwise, try to add it to the bucket.
					if (iBucket.isEmpty())
						this.farCount++;
					else {
						this.storeGenome(seq, iBucket.getAsInt());
					}
					inCount++;
					// Check for a progress message.
					long now = System.currentTimeMillis();
					if (now - lastMsg >= 10000) {
						log.info("{} genomes processed, {} kept, {} too far, {} full buckets, {} empty buckets.",
								inCount, this.keptCount, this.farCount, this.fullCount, this.emptyCount);
						lastMsg = now;
					}
				}
			}
			log.info("All done. {} genomes processed, {} kept, {} too far, {} full buckets, {} empty buckets.",
					inCount, this.keptCount, this.farCount, this.fullCount, this.emptyCount);
		}
		// Now we write the output. We need to combine all the buckets and sort them.
		Map<String, String> allMap = new TreeMap<String, String>();
		this.bucketList.stream().forEach(x -> allMap.putAll(x));
		log.info("Writing {} genomes to output.", allMap.size());
		writer.println("genome_id\tspecies");
		for (var gEntry : allMap.entrySet())
			writer.println(gEntry.getKey() + "\t" + gEntry.getValue());
	}

	/**
	 * This method reads in the first sequence, saves its kmer object, and stores it in the first
	 * bucket.
	 *
	 * @param inStream	FASTA input stream from the finder file
	 */
	private void readFirstSequence(FastaInputStream inStream) {
		Sequence firstSeq = inStream.next();
		this.sourceKmers = new DnaKmers(firstSeq.getSequence());
		// Store the sequence in the first bucket.
		String genome1 = this.storeGenome(firstSeq, 0);
		log.info("First sequence processed from genome {}.", genome1);
		// Save the primer ID.
		this.primerID = firstSeq.getLabel();
	}

	/**
	 * Store a genome in a bucket. The sequence contains a feature ID from which the genome ID can
	 * be computed as the label. The comment contains the species name after the first tab. (This
	 * is a characteristic of finder FASTA files. The DNA itself is not needed.
	 *
	 * @param seq		genome sequence to store.
	 * @param idx		index of the target bucket
	 *
	 * @return the ID of the genome added
	 */
	private String storeGenome(Sequence seq, int idx) {
		// Get the genome ID and species.
		String retVal = Feature.genomeOf(seq.getLabel());
		String genomeName = StringUtils.substringAfter(seq.getComment(), "\t");
		// Figure out if we have a bucket.
		Map<String, String> bucket = this.bucketList.get(idx);
		if (bucket.isEmpty()) {
			this.emptyCount--;
			this.keptCount++;
			bucket.put(retVal, genomeName);
		} else {
			// Get the current bucket size.
			int oldSize = bucket.size();
			// Only proceed if the bucket has room.
			if (oldSize < this.bucketSize) {
				// Add this genome. We have a lot of bookkeeping to do. We need to know
				// if the genome is new or duplicate (duplicates are rare, but happen).
				// We also need to know if the bucket is now full.
				bucket.put(retVal, genomeName);
				int newSize = bucket.size();
				if (newSize > oldSize) {
					this.keptCount++;
					if (newSize >= this.bucketSize)
						this.fullCount++;
				}
			}
		}
		return retVal;
	}

}
