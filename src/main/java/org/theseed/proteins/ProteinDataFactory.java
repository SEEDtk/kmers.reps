/**
 *
 */
package org.theseed.proteins;

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import java.util.Iterator;

import org.theseed.genome.Feature;
import org.theseed.p3api.Connection;
import org.theseed.p3api.Connection.Table;

import com.github.cliftonlabs.json_simple.JsonObject;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class is responsible for creating and managing ProteinData objects.  It
 * maintains a genus/species table and a connection to PATRIC.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinDataFactory implements Iterable<ProteinData> {

    public static final String SEED_FUNCTION = "Phenylalanyl-tRNA synthetase alpha chain";

    // FIELDS

    /** logging facility */
    private static Logger log = LoggerFactory.getLogger(ProteinDataFactory.class);
    /** set of taxon IDs associated with genera */
    private Set<String> genera;
    /* map of taxon IDs associated with species to genetic codes */
    private HashMap<String, Integer> species;
    /** connection to PATRIC */
    private Connection p3;
    /** master list of ProteinData objects */
    List<ProteinData> master;
    /** map of ProteinData objects by genome ID */
    Map<String, ProteinData> idMap;

    /**
     * Initialize the PATRIC connection and read in the genus and species sets.
     */
    public ProteinDataFactory() {
        // Connect to PATRIC.
        this.p3 = new Connection();
        // Get all the species.
        log.info("Downloading taxonomic data.");
        List<JsonObject> taxonList = p3.query(Table.TAXONOMY,
                "taxon_id,genetic_code", "eq(taxon_rank,species)");
        this.species = new HashMap<String, Integer>(taxonList.size());
        for (JsonObject taxonData : taxonList) {
            String speciesId = Connection.getString(taxonData, "taxon_id");
            int gc = Connection.getInt(taxonData, "genetic_code");
            this.species.put(speciesId, gc);
        }
        log.info("{} species tabulated.", this.species.size());
        // Get all the genera.
        taxonList = p3.query(Table.TAXONOMY, "taxon_id", "eq(taxon_rank,genus)");
        this.genera = taxonList.stream().map(x -> Connection.getString(x, "taxon_id")).collect(Collectors.toSet());
        log.info("{} genera tabulated.", this.genera.size());
        // Create the master list and map.
        this.master = new ArrayList<ProteinData>(200000);
        this.idMap = new HashMap<String, ProteinData>(200000);
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
     * @return TRUE if successful, FALSE if the genus or species could not be determined
     */
    public boolean addGenome(String genomeId, String genomeName, String lineageString, double score) {
        // Parse out the lineage.
        String genus = null;
        String species = null;
        String domain = "Bacteria";
        int gc = 11;
        String[] lineage = StringUtils.splitByWholeSeparator(lineageString, "::");
        for (String taxId : lineage) {
            if (this.genera.contains(taxId))
                genus = taxId;
            else if (taxId.contentEquals("2157"))
                domain = "Archaea";
            else {
                Integer gCode = this.species.get(taxId);
                if (gCode != null) {
                    gc = gCode;
                    species = taxId;
                }
            }
        }
        boolean retVal = (genus != null && species != null);
        // If we found the genus and species, create the protein data object.
        if (retVal) {
            ProteinData newGenome = new ProteinData(genomeId, genomeName, domain, genus, species,
                    gc, score);
            this.master.add(newGenome);
            this.idMap.put(genomeId, newGenome);
        } else {
            log.debug("Missing genus or species for {}.", genomeId);
        }
        return retVal;
    }

    /**
     * Finish the master list.  The list is traversed sequentially, retrieving the DNA and protein
     * sequences.  This is a slow, expensive operation requiring a lot of PATRIC access.  The
     * genus and species sets are released at the beginning, since they are no longer needed and
     * we are about to use a lot of memory for the sequences.
     *
     * @param batchSize		number of genomes to process in each batch
     *
     * @throws UnsupportedEncodingException
     */
    public void finishList(int batchSize) throws UnsupportedEncodingException {
        this.genera = null;
        this.species = null;
        // We ask for all the features relating to the seed protein function.  We use a special
        // role map to isolate it.
        RoleMap seedMap = new RoleMap();
        seedMap.findOrInsert(SEED_FUNCTION);
        // Ask for all the features with this function.
        log.info("Retrieving seed protein features.");
        List<JsonObject> features = p3.query(Table.FEATURE,
                "genome_id,patric_id,product,aa_sequence_md5,na_sequence_md5",
                "eq(product," + URLEncoder.encode(SEED_FUNCTION, StandardCharsets.UTF_8.toString()) + ")",
                "eq(annotation,PATRIC)");
        // These maps are keyed by MD5, and map each MD5 to the list of ProteinData objects for the
        // associated genomes.
        Map<String, Collection<ProteinData>> dnaMap = new HashMap<String, Collection<ProteinData>>(batchSize);
        Map<String, Collection<ProteinData>> protMap = new HashMap<String, Collection<ProteinData>>(batchSize);
        // We are ready.  Loop through the features, retrieving the sequences.
        log.info("Retrieving DNA and protein sequences.");
        for (JsonObject feature : features) {
            // Check this feature for a valid function.  It usually IS valid.  Rarely, we get a substring match of something that
            // is similar, but not correct.
            String[] roleNames = Feature.rolesOfFunction(Connection.getString(feature, "product"));
            boolean foundRole = Arrays.stream(roleNames).anyMatch(x -> seedMap.containsName(x));
            if (foundRole) {
                // Get the protein data for the feature's genome.
                String genomeId = Connection.getString(feature, "genome_id");
                ProteinData genomeData = this.idMap.get(genomeId);
                // Only proceed if we found it.  If we didn't find it, then it is not one of our genomes.
                if (genomeData != null) {
                    // Verify that we have a valid feature ID and both MD5s.  Note that there is no trace message for a missing
                    // feature ID, as features with a missing ID have a special meaning.
                    String fid = Connection.getString(feature, "patric_id");
                    String dnaMd5 = Connection.getString(feature, "na_sequence_md5");
                    String protMd5 = Connection.getString(feature, "aa_sequence_md5");
                    if (dnaMd5 == null || dnaMd5.isEmpty()) {
                        log.debug("Missing DNA sequence for seed protein of {}.", genomeId);
                    } else if (protMd5 == null || protMd5.isEmpty()) {
                        log.debug("Missing protein sequence for seed protein of {}.", genomeId);
                    } else if (fid != null && ! fid.isEmpty()) {
                        genomeData.setFid(fid);
                        dnaMap.computeIfAbsent(dnaMd5, k -> new ArrayList<ProteinData>(5)).add(genomeData);
                        protMap.computeIfAbsent(protMd5, k -> new ArrayList<ProteinData>(5)).add(genomeData);
                        // If this fills a batch, process it.
                        if (dnaMap.size() >= batchSize) {
                            this.processMaps(dnaMap, protMap);
                        }
                    }
                }
            }
        }
        // Process the residual batch.
        if (dnaMap.size() > 0) this.processMaps(dnaMap, protMap);
        // Now run through and remove the genomes that aren't filled in.  Move the
        // genomes with ambiguity characters to the end.
        List<ProteinData> deferList = new ArrayList<ProteinData>(1000);
        log.info("Removing genomes with incomplete data and resorting list.");
        int deleteCount = 0;
        Iterator<ProteinData> iter = this.master.iterator();
        while (iter.hasNext()) {
            ProteinData genomeData = iter.next();
            if (genomeData.getDna() == null || genomeData.getProtein() == null) {
                iter.remove();
                this.idMap.remove(genomeData.getGenomeId());
                deleteCount++;
            } else if (genomeData.getProtein().contains("X")) {
                iter.remove();
                deferList.add(genomeData);
            }
        }
        log.info("{} incomplete genomes removed. {} ambiguous proteins relocated.",
                deleteCount, deferList.size());
        // Add the deferred proteins to the end of the list.
        this.master.addAll(deferList);
    }

    /**
     * Here the accumulated protein and DNA MD5 maps are used to query PATRIC for the actual protein and DNA sequences.
     *
     * @param dnaMap	map of DNA MD5s to ProteinData objects
     * @param protMap	map of protein MD%s to ProteinData objects
     */
    private void processMaps(Map<String, Collection<ProteinData>> dnaMap, Map<String, Collection<ProteinData>> protMap) {
        log.info("Retrieving DNA sequences for {} proteins.", dnaMap.size());
        Map<String, JsonObject> sequences = p3.getRecords(Table.SEQUENCE, dnaMap.keySet(), "sequence");
        int dnaSet = 0;
        int protSet = 0;
        for (Map.Entry<String, JsonObject> sequence : sequences.entrySet()) {
            Collection<ProteinData> genomeData = dnaMap.get(sequence.getKey());
            String dna = Connection.getString(sequence.getValue(), "sequence");
            for (ProteinData genomeDatum : genomeData) {
                genomeDatum.setDna(dna);
                dnaSet++;
            }
        }
        log.info("Retrieving protein sequences for {} proteins.", protMap.size());
        sequences = p3.getRecords(Table.SEQUENCE, protMap.keySet(), "sequence");
        for (Map.Entry<String, JsonObject> sequence : sequences.entrySet()) {
            Collection<ProteinData> genomeData = protMap.get(sequence.getKey());
            String prot = Connection.getString(sequence.getValue(), "sequence");
            for (ProteinData genomeDatum : genomeData) {
                genomeDatum.setProtein(prot);
                protSet++;
            }
        }
        log.info("{} DNA sequences and {} protein sequences stored.", dnaSet, protSet);
        // Erase the maps so they can be refilled for the next batch.
        dnaMap.clear();
        protMap.clear();
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

}
