/**
 *
 */
package org.theseed.proteins.kmers;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.p3api.P3TaxData;
import org.theseed.p3api.SolrFilter;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.roles.RoleUtilities;
import org.theseed.sequence.DnaKmers;
import org.theseed.stats.GenomeEval;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This class is responsible for creating and managing ProteinData objects.  It
 * maintains a genus/species table and a connection to PATRIC.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinDataFactory implements Iterable<ProteinData> {

    // FIELDS

    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(ProteinDataFactory.class);
    /** taxonomy map */
    private P3TaxData taxMap;
    /** connection to PATRIC */
    private final P3CursorConnection p3;
    /** master list of ProteinData objects */
    private final SortedSet<ProteinData> master;
    /** map of ProteinData objects by genome ID */
    private final Map<String, ProteinData> idMap;
    /** map of special NCBI genomes to ratings */
    protected Map<String, ProteinData.Rating> ncbiRefMap;
    /** genome statistics */
    private final GenomeEval statistics;
    /** list of supported repgen levels */
    private int[] repLevels;
    /** minimum SSU rRNA length */
    public static final int MIN_SSU_LEN = 1400;
    /** minimum SSU rRNA length for validation */
    public static final int USEFUL_SSU_LEN = 700;
    /** maximum SSU rRNA distance */
    public static final double MAX_SSU_DISTANCE = 0.5;
    /** maximum SSU rRNA distance within a genus*/
    public static final double MAX_GENUS_SSU_DISTANCE = 0.75;
    /** seed protein function */
    public static final String SEED_FUNCTION = "Phenylalanyl-tRNA synthetase alpha chain";

    /**
     * Initialize the PATRIC connection and read in the genus and species sets.
     */
    public ProteinDataFactory() {
        // Connect to PATRIC.
        this.p3 = new P3CursorConnection();
        // Get all the taxonomy stuff.
        this.taxMap = new P3TaxData(this.p3);
        // Get the NCBI reference and representative genomes.
        List<JsonObject> ncbiGenomes;
        try {
            ncbiGenomes = p3.getRecords("genome", P3CursorConnection.MAX_LIMIT, "genome_id,reference_genome",
                    SolrFilter.IN("reference_genome", "Representative", "Reference"));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        log.info("{} NCBI special genomes found.", ncbiGenomes.size());
        this.ncbiRefMap = new HashMap<>(ncbiGenomes.size());
        for (JsonObject genome : ncbiGenomes) {
            String genomeId = KeyBuffer.getString(genome, "genome_id");
            if (genomeId != null && ! genomeId.isEmpty()) {
                String type = KeyBuffer.getString(genome, "reference_genome");
                if (type != null) {
                    if (type.contentEquals("Reference"))
                        this.ncbiRefMap.put(genomeId, ProteinData.Rating.NCBI_REF);
                    else if (type.contentEquals("Representative"))
                        this.ncbiRefMap.put(genomeId, ProteinData.Rating.NCBI_REP);
                }
            }
        }
        // Create the master list and map.
        this.master = new TreeSet<>();
        this.idMap = new HashMap<>(200000);
        // Create the statistics object.
        this.statistics = new GenomeEval();
        // Denote we have no repgen levels yet.
        this.repLevels = new int[0];
    }

    /**
     * Add a genome to the end of the master list.
     *
     * @param genomeId			ID of the new genome
     * @param genomeName		name of the new genome
     * @param lineageString		taxonomic lineage string for the new genome (IDs separated by
     * 							double colons)
     * @param score				quality score
     *
     * @return TRUE if successful, FALSE if the species could not be determined
     */
    public boolean addGenome(String genomeId, String genomeName, String lineageString, double score) {
        // Parse out the lineage.
        String genus = null;
        String species = null;
        String domain = "Bacteria";
        int gc = 11;
        String[] lineage = StringUtils.splitByWholeSeparator(lineageString, "::");
        for (String taxId : lineage) {
            if (this.taxMap.isGenus(taxId))
                genus = taxId;
            else if (taxId.contentEquals("2157"))
                domain = "Archaea";
            else {
                int gCode = this.taxMap.checkSpecies(taxId);
                if (gCode != 0) {
                    gc = gCode;
                    species = taxId;
                }
            }
        }
        boolean retVal = (species != null);
        // If we found the species, create the protein data object.
        if (retVal) {
            ProteinData newGenome = new ProteinData(genomeId, genomeName, domain, genus, species,
                    gc, score);
            // Initialize the rating to one of the three good ones.
            newGenome.setRating(this.ncbiRefMap.getOrDefault(genomeId, ProteinData.Rating.NORMAL));
            // Put the genome in the data structures.
            this.master.add(newGenome);
            this.idMap.put(genomeId, newGenome);
        } else {
            // Otherwise, note that the species is missing.
            log.debug("Missing species for {}.", genomeId);
            this.statistics.count("MissingSpecies", 1);
        }
        return retVal;
    }

    /**
     * Analyze a summary file line to update the statistics.
     *
     * @param line	summary file line to analyze
     */
    public void analyze(TabbedLineReader.Line line) {
        this.statistics.analyze(line);
    }

    /**
     * Finish the master list.  The list is traversed twice sequentially, retrieving the DNA and protein
     * sequences.  This is a slow, expensive operation requiring a lot of PATRIC access.  The
     * genus and species sets are released at the beginning, since they are no longer needed and
     * we are about to use a lot of memory for the sequences.
     *
     * Note that for the seed protein, we will only find one, because a good genome has only one.
     * For the SSU rRNAs, we may find many, so we need a different strategy.
     *
     * @param batchSize		number of genomes to process in each batch
     *
     * @throws UnsupportedEncodingException
     */
    public void finishList(int batchSize) throws UnsupportedEncodingException {
        this.taxMap = null;
        // We ask for all the features relating to the seed protein function.  We use a special
        // role map to isolate it.
        RoleMap seedMap = new RoleMap();
        seedMap.findOrInsert(SEED_FUNCTION);
        // Ask for all the features with this function in the specified genomes.
        log.info("Retrieving seed protein features.");
        this.getFeaturesByFunction(batchSize, seedMap);
        // Now run through and remove the genomes that aren't filled in or have multiple
        // ambiguity characters in the seed protein.
        log.info("Removing genomes with incomplete data or ambiguity.");
        int deleteCount = 0;
        Iterator<ProteinData> iter = this.master.iterator();
        while (iter.hasNext()) {
            ProteinData genomeData = iter.next();
            boolean removeFlag = false;
            if (genomeData.getDna() == null) {
                removeFlag = true;
                this.statistics.count("Missing Seed DNA", 1);
            }
            if (genomeData.getProtein() == null) {
                removeFlag = true;
                this.statistics.count("Missing Seed AA", 1);
            } else if (genomeData.getProtein().contains("XX")) {
                removeFlag = true;
                this.statistics.count("Ambiguous Seed", 1);
            }
            if (removeFlag) {
                iter.remove();
                this.idMap.remove(genomeData.getGenomeId());
                deleteCount++;
            }
        }
        log.info("{} incomplete or ambiguous genomes removed.",
                deleteCount);
        // Now we do the SSU rRNAs with our slightly-reduced genome set.
        Map<String, List<JsonObject>> genomeMap = new HashMap<>(this.master.size());
        this.retrieveSsuRnaFeatures(genomeMap);
        // We have all the SSU rRNA sequences, which is quite an accomplishment.  Run through the
        // genomes, collecting them and updating the SSUs.  This involves changing the key to a
        // sorted set, so we clear the master and rebuild it.
        List<ProteinData> masterClone = new ArrayList<>(this.master);
        this.master.clear();
        // For single-RNA genomes, this tracks good-RNA genomes for each genus.
        Map<String, DnaKmers> genusMap = new HashMap<>(10000);
        List<ProteinData> retryQueue = new ArrayList<>(masterClone.size());
        // We need some counters to display.
        deleteCount = 0;
        int badCount = 0;
        int shortCount = 0;
        int processCount = 0;
        for (ProteinData genomeData : masterClone) {
            String genomeId = genomeData.getGenomeId();
            List<JsonObject> jsons = genomeMap.get(genomeId);
            if (jsons == null) {
                // Here we can't find a 16s.  Reject the genome.
                log.warn("No 16S sequence found for {}: {}.", genomeId, genomeData.getGenomeName());
                deleteCount++;
            } else {
                // Get all the sequences for this genome.
                List<String> rnas = jsons.stream().map(x -> KeyBuffer.getString(x, "na_sequence"))
                        .filter(x -> x != null && ! x.isEmpty()).collect(Collectors.toList());
                if (rnas.isEmpty()) {
                    // Here the SSUs had no valid sequences.  Reject the genome.
                    log.warn("All 16s sequences invalid in {}: {}.", genomeId, genomeData.getGenomeName());
                    deleteCount++;
                    this.statistics.count("16s Seq Not Found", 1);
                } else {
                    // We have sequences.  Merge them into the protein data.  This may change the rating.
                    ProteinData.Rating rating = genomeData.setSsuSequence(rnas);
                    String genus = genomeData.getGenus();
                    switch (rating) {
                    case BAD_SSU -> {
                        badCount++;
                        this.master.add(genomeData);
                        }
                    case SHORT_SSU -> {
                        shortCount++;
                        this.master.add(genomeData);
                        }
                    case SINGLE_SSU -> {
                        // Here we want to queue the genome for a second pass if it has a genus.
                        if (genus != null)
                            retryQueue.add(genomeData);
                        }
                    default -> {
                        // Here the genome is good.  If it is the first for this genus, save it
                        // for the retry queue.
                        if (genus != null)
                            genusMap.putIfAbsent(genus, new DnaKmers(genomeData.getSsuSequence()));
                        this.master.add(genomeData);
                        }
                    }
                }
            }
            processCount++;
            if (log.isInfoEnabled() && processCount % 2000 == 0)
                log.info("{} genomes processed, {} bad, {} short, {} queued for retry.", processCount,
                        badCount, shortCount, retryQueue.size());
        }
        log.info("{} genomes deleted due to missing sequences.  {} flagged as bad, {} as short",
                deleteCount, badCount, shortCount);
        // Now we try once again for the single-SSU genomes.  We are hoping that we can use the reference
        // genome for the genus to verify the single SSU.
        log.info("Retrying {} genomes at genus level.", retryQueue.size());
        int reclaimed = 0;
        int processed = 0;
        for (ProteinData genomeData : retryQueue) {
            DnaKmers genusKmers = genusMap.get(genomeData.getGenus());
            // The rating at this point is SINGLE_SSU which is worse than short.
            if (genusKmers != null) {
                // Here there is data we can use to reclaim the genome.
                DnaKmers genomeKmers = new DnaKmers(genomeData.getSsuSequence());
                double dist = genusKmers.distance(genomeKmers);
                ProteinData.Rating rating;
                if (dist >= MAX_GENUS_SSU_DISTANCE) {
                    // The SSU is outside the genus, so it is bad.
                    rating = ProteinData.Rating.BAD_SSU;
                } else {
                    // Here we can reclaim it, but it might still be short.
                    if (genusKmers.getLen() < MIN_SSU_LEN)
                        rating = ProteinData.Rating.SHORT_SSU;
                    else
                        rating = ProteinData.Rating.NORMAL;
                    reclaimed++;
                }
                genomeData.setRating(rating);
            }
            this.master.add(genomeData);
            processed++;
            if (log.isInfoEnabled() && this.master.size() % 5000 == 0)
                log.info("{} genomes checked.  {} reclaimed.", processed, reclaimed);
        }
        log.info("Final genome count is {}.  {} were reclaimed by genus testing.",
                this.master.size(), reclaimed);
    }

    /**
     * Retrieve SSU rRNA features from the PATRIC database into the specified maps.
     * Form the features into a map based on genome ID.  Also, we create a map that will
     * list all the SSU rRNA features for each genome.
     *
     * @param genomeMap     a map from genome ID to list of SSU rRNA feature records
     */
    private void retrieveSsuRnaFeatures(Map<String, List<JsonObject>> genomeMap) {
        log.info("Retrieving SSU rRNA features for {} genomes.", this.idMap.size());
        List<JsonObject> features;
        try {
            features = p3.getRecords("feature", P3CursorConnection.MAX_LIMIT, 2000, "genome_id", this.idMap.keySet(),
                    "genome_id,patric_id,product,na_sequence",
                    SolrFilter.IN("feature_type", "rrna", "rRNA", "misc_RNA"),
                    SolrFilter.EQ("annotation", "PATRIC"));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        log.info("{} total rRNAs found.", features.size());
        int found = 0;
        for (JsonObject feature : features) {
            String product = KeyBuffer.getString(feature, "product");
            if (RoleUtilities.SSU_R_RNA.matcher(product).find()) {
                // Here we have a real 16S feature.
                String genomeId = KeyBuffer.getString(feature, "genome_id");
                List<JsonObject> feats = genomeMap.computeIfAbsent(genomeId, x -> new ArrayList<JsonObject>(5));
                feats.add(feature);
                found++;
            }
        }
        log.info("{} 16s features found in {} genomes.", found, genomeMap.size());
    }

    /**
     * Locate all of the features with the seed function.
     *
     * @param batchSize     batch size for PATRIC queries
     * @param seedMap       a role map containing the seed functions
     */
    private void getFeaturesByFunction(int batchSize, RoleMap seedMap) {
        log.info("Querying PATRIC for seed protein features in {} genomes.", this.idMap.size());
        List<JsonObject> features;
        try {
            features = p3.getRecords("feature", P3CursorConnection.MAX_LIMIT, batchSize, "genome_id", this.idMap.keySet(),
                    "genome_id,patric_id,product,na_sequence,aa_sequence",
                    SolrFilter.EQ("product", SEED_FUNCTION));
        } catch (IOException e) {
            throw new UncheckedIOException(e);       }
        // We are ready.  Loop through the features, retrieving the sequences.
        log.info("{} seed protein features found.", features.size());
        int badRole = 0;
        int badDna = 0;
        int badProt = 0;
        int badFid = 0;
        for (JsonObject feature : features) {
            // Check this feature for a valid function.  It usually IS valid.  Rarely, we get a substring match of something that
            // is similar, but not correct.
            String[] roleNames = Feature.rolesOfFunction(KeyBuffer.getString(feature, "product"));
            boolean foundRole = Arrays.stream(roleNames).anyMatch(x -> seedMap.containsName(x));
            if (! foundRole) {
                log.debug("Seed protein feature {} has invalid function {}.", KeyBuffer.getString(feature, "patric_id"),
                        KeyBuffer.getString(feature, "product"));
                badRole++;
            } else {
                // Get the protein data for the feature's genome.
                String genomeId = KeyBuffer.getString(feature, "genome_id");
                ProteinData genomeData = this.idMap.get(genomeId);
                // Only proceed if we found it.  If we didn't find it, then it is not one of our genomes.
                if (genomeData != null) {
                    // Verify that we have a valid feature ID and both sequences.  Note that there is no trace message for a missing
                    // feature ID, as features with a missing ID have a special meaning.
                    String fid = KeyBuffer.getString(feature, "patric_id");
                    String dna = KeyBuffer.getString(feature, "na_sequence");
                    String prot = KeyBuffer.getString(feature, "aa_sequence");
                    if (dna == null || dna.isEmpty()) {
                        log.debug("Missing DNA sequence for seed protein of {}.", genomeId);
                        badDna++;
                    } else if (prot == null || prot.isEmpty()) {
                        log.debug("Missing protein sequence for seed protein of {}.", genomeId);
                        badProt++;
                    } else if (fid == null || ! Feature.FID_PARSER.matcher(fid).matches()) {
                        log.debug("Missing or invalid feature ID for seed protein of {}.", genomeId);
                        badFid++;
                    } else {
                        genomeData.setFid(fid);
                        genomeData.setDna(dna);
                        genomeData.setProtein(prot);
                    }
                }
            }
        }
        log.info("{} invalid functions, {} missing DNA, {} missing protein, {} missing feature ID.", badRole, badDna, badProt, badFid);
    }

    @Override
    public Iterator<ProteinData> iterator() {
        return this.master.iterator();
    }

    /**
     * @return the genome with a specific ID
     */
    public ProteinData getGenome(String genomeId) {
        return this.idMap.get(genomeId);
    }

    /**
     * @return the number of genomes in this list
     */
    public int size() {
        return this.master.size();
    }

    /**
     * Restore the repgen sets from the input directory.  This insures that the old
     * representatives are kept when we rebuild everything.
     *
     * @param container		container in which to store the repgen sets
     * @param inDir			input directory
     *
     * @throws IOException
     */
    public static void restoreData(IRepGenContainer container, File inDir) throws IOException {
        // Find the repXX.ser files.
        File[] repDbFiles = inDir.listFiles((d, n) -> n.matches("rep\\d+\\.ser"));
        // Initialize the repgen set list.
        container.initRepGenSets(repDbFiles.length);
        // Loop through the files.
        for (File repDbFile : repDbFiles) {
            log.info("Loading repgen set from {}.", repDbFile);
            RepGenomeDb repDb = RepGenomeDb.load(repDbFile);
            container.addRepGenSet(repDb);
        }
    }

    /**
     * Remove genomes from the list with a low rating.
     *
     * @param rating	minimum rating to keep
     */
    public void prune(ProteinData.Rating rating) {
        // We create a tiny little map for the rating names.
        Map<ProteinData.Rating, String> ratingNames = new EnumMap<>(ProteinData.Rating.class);
        Arrays.stream(ProteinData.Rating.values())
                .forEach(x -> ratingNames.put(x, "Rating-" + x.toString()));
        int deleteCount = 0;
        Iterator<ProteinData> iter = this.master.iterator();
        while (iter.hasNext()) {
            ProteinData curr = iter.next();
            ProteinData.Rating currR = curr.getRating();
            // Count the genome.
            this.statistics.count(ratingNames.get(currR), 1);
            // Remove the genome if it is rated too low.
            if (curr.getRating().compareTo(rating) > 0) {
                iter.remove();
                deleteCount++;
            }
        }
        log.info("{} genomes with rating less than {} deleted.", deleteCount, rating);
        this.statistics.count("Low Rating", deleteCount);
        this.statistics.count("Final Genome Count", this.master.size());
    }

    /**
     * Write the statistics to the specified file.
     *
     * @param outfile	output file name
     *
     * @throws IOException
     */
    public void writeStats(File outFile) throws IOException {
        this.statistics.write(outFile);
    }

    /**
     * @return TRUE if the specified genome is in this factory
     *
     * @param genomeId		ID of the genome to check
     */
    public boolean contains(String genomeId) {
        return this.idMap.containsKey(genomeId);
    }

    /**
     * @return the array of supported repgen levels
     */
    public int[] getRepLevels() {
        return this.repLevels;
    }

    /**
     * Store the repgen set levels.
     *
     * @param levels	array of repgen levels to store
     */
    public void setRepLevels(int[] levels) {
        this.repLevels = Arrays.copyOf(levels, levels.length);
    }

    /**
     * Write a report on the genomes in the protein data to the specified file.  This is our
     * best list of genome quality data, since it includes the SSU analysis.
     *
     * @param file	output file for the report
     *
     * @throws IOException
     */
    public void writeReport(File file) throws IOException {
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(ProteinData.getHeader(this));
            for (var protein : this.master)
                writer.println(protein.getLine(this));
        }
    }

}
