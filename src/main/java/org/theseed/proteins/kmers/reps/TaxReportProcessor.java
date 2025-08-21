/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;

/**
 * This command processes the good-genome list to determine what a particular taxonomic ID tells us about membership in a repgen set.
 * In particular, for each genus and species, we compute what fraction of each repgen set is in that taxonomic group.  To do this,
 * we associate the taxon IDs with genomes from the good-genome list (which includes the genus and species IDs for each genome as
 * well as its representative genome ID in each repgen set).  The weight of the taxon ID for a repgen group is the number
 * of genomes in that taxon belonging to the group divided by the total number of genomes in that taxon.  The output will
 * be a table of taxon ID / repgen ID pairs listing the weights.
 *
 * The standard input should contain the appropriate repgen list file.  The report will be written to the standard output.
 * The positional parameter is the name of the appropriate repgen group ("repXX"), which will be a column heading in the
 * good-genome list.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input good-genome list file (if not STDIN)
 * -o	output report file (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class TaxReportProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(TaxReportProcessor.class);
    /** map of taxon IDs to count maps, each count map keyed by repgen ID */
    private Map<Integer, CountMap<String>> taxCountMaps;
    /** map of taxon IDs to ranks (genus or species) */
    private Map<Integer, String> rankMap;
    /** index of the genus ID input column */
    private int genusColIdx;
    /** index of the species ID input column */
    private int speciesColIdx;
    /** index of the repgen ID column */
    private int repColIdx;
    /** number of missing taxonomic IDs */
    private int missingCount;

    // COMMAND-LINE OPTIONS

    /** repgen set name */
    @Argument(index = 0, metaVar = "repXX", usage = "name of the input column containing the repgen IDs of interest", required = true)
    private String repColId;

    @Override
    protected void setPipeDefaults() {
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        log.info("Analyzing column headings for input.");
        // Get the genus and species columns.
        this.genusColIdx = inputStream.findField("genus");
        this.speciesColIdx = inputStream.findField("species");
        // Find the representative-genome ID column.
        this.repColIdx = inputStream.findField(repColId);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Create the master count map hash and the rank hash.
        this.taxCountMaps = new HashMap<Integer, CountMap<String>>(10000);
        this.rankMap = new HashMap<Integer, String>(10000);
        // Clear the missing-ID counter.
        this.missingCount = 0;
        // Count the genomes processed.
        int gCount = 0;
        // We loop through the input file, counting genomes.
        log.info("Processing input file.");
        for (var line : inputStream) {
            int genus = line.getInt(genusColIdx);
            int species = line.getInt(speciesColIdx);
            // Get the repgen group ID.
            String repgenId = line.get(repColIdx);
            // Store the genus and species.
            this.storeTaxon(genus, repgenId, "genus");
            this.storeTaxon(species, repgenId, "species");
            // Count this genome.
            gCount++;
            if (log.isInfoEnabled() && gCount % 10000 == 0)
                log.info("{} genomes processed. {} taxonomic groupings found.", gCount, this.rankMap.size());
        }
        log.info("{} total genomes processed.  {} taxonomic IDs were missing.", gCount, this.missingCount);
        // Now we want to write the output report.  For each taxonomic grouping, we output the weight for each repgen ID.
        writer.println("tax_id\trank\trepgen_id\tcount\tweight");
        for (var taxCountMap : this.taxCountMaps.entrySet()) {
            CountMap<String> repCounts = taxCountMap.getValue();
            int taxId = taxCountMap.getKey();
            String taxIdString = Integer.toString(taxId);
            String rank = this.rankMap.getOrDefault(taxId, "<unknown>");
            // Get the sum of the counts for this group.
            double sum = repCounts.getTotal();
            var countList = repCounts.sortedCounts();
            for (var countData : countList) {
                int count = countData.getCount();
                double weight = count / sum;
                String repgenId = countData.getKey();
                writer.println(taxIdString + "\t" + rank + "\t" + repgenId + "\t" + Integer.toString(count)
                        + "\t" + Double.toString(weight));
            }
        }
    }

    /**
     * Count a represented genome in a specific taxon.
     *
     * @param taxId		ID of the taxonomic grouping containing the genome
     * @param repgenId	ID of the genome's representative
     * @param rank		rank of the taxonomic grouping ("genus" or "species")
     */
    private void storeTaxon(int taxId, String repgenId, String rank) {
        // Check for a missing value.
        if (taxId == 0)
            this.missingCount++;
        else {
            this.rankMap.put(taxId, rank);
            CountMap<String> repCounts = this.taxCountMaps.computeIfAbsent(taxId, x -> new CountMap<String>());
            repCounts.count(repgenId);
        }
    }

}
