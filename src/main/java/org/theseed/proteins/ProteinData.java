/**
 *
 */
package org.theseed.proteins;

import java.util.HashMap;
import java.util.Map;

import org.theseed.proteins.kmers.ProteinKmers;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.proteins.kmers.reps.RepGenomeDb.Representation;

/**
 * This class describes the seed protein for a genome.  It contains the DNA and amino acid sequences.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinData {

    /** feature ID */
    private String fid;
    /** genome ID */
    private String genomeId;
    /** genome name */
    private String genomeName;
    /** DNA sequence */
    private String dna;
    /** amino acid sequence */
    private String protein;
    /** taxonomic genus */
    private String genus;
    /** taxonomic species */
    private String species;
    /** quality score */
    private double score;
    /** domain name */
    private String domain;
    /** genetic code */
    private int geneticCode;
    /** proteinkmers object for computing representation */
    private ProteinKmers kmers;
    /** representative genome information for each score level */
    private Map<Integer, RepGenomeDb.Representation> representation;

    /**
     * Create a new protein data object.
     *
     * @param genomeId		ID of the genome
     * @param genomeName	name of the genome
     * @param domain		name of the domain-- Bacteria or Archaea
     * @param genus			taxonomic ID of the genus
     * @param species		taxonomic ID of the species
     * @param gc			genetic code
     * @param score			quality score
     */
    public ProteinData(String genomeId, String genomeName, String domain, String genus, String species,
            int gc, double score) {
        this.genomeId = genomeId;
        this.genomeName = genomeName;
        this.genus = genus;
        this.species = species;
        this.geneticCode = gc;
        this.score = score;
        this.domain = domain;
        this.representation = new HashMap<Integer, RepGenomeDb.Representation>();
        this.kmers = null;
    }
    /**
     * @return the ID of the seed protein feature
     */
    public String getFid() {
        return fid;
    }
    /**
     * @return the ID of this genome
     */
    public String getGenomeId() {
        return genomeId;
    }
    /**
     * @return the name of this genome
     */
    public String getGenomeName() {
        return genomeName;
    }
    /**
     * @return the dna sequence of the seed protein
     */
    public String getDna() {
        return dna;
    }
    /**
     * @return the protein
     */
    public String getProtein() {
        return protein;
    }
    /**
     * @return the genus
     */
    public String getGenus() {
        return genus;
    }
    /**
     * @return the species
     */
    public String getSpecies() {
        return species;
    }

    /**
     * Store the seed protein feature ID.
     *
     * @param fid	feature ID to store
     */
    public void setFid(String fid) {
        this.fid = fid;
    }
    /**
     * @return the quality score
     */
    public double getScore() {
        return score;
    }
    /**
     * Store the seed protein DNA string.
     *
     * @param dna 	the dna string to store
     */
    public void setDna(String dna) {
        this.dna = dna;
    }
    /**
     * Store the seed protein amino acid string.
     *
     * @param protein 	the protein string to store
     */
    public void setProtein(String protein) {
        this.protein = protein;
    }
    /**
     * @return the domain
     */
    public String getDomain() {
        return domain;
    }
    /**
     * @return the genetic code of this genome
     */
    public int getGeneticCode() {
        return geneticCode;
    }

    /**
     * @return TRUE if this genome has a representative for the specified repgen threshold
     *
     * @param sim	similarity threshold to check
     */
    public boolean checkRepresented(int sim) {
        return this.representation.containsKey(sim);
    }
    public Representation getRepresentation(RepGenomeDb repGenSet) {
        int sim = repGenSet.getThreshold();
        RepGenomeDb.Representation retVal = this.representation.get(sim);
        // If we don't already have a representative, compute one.
        if (retVal == null) {
            // Insure we have the protein kmers.
            if (this.kmers == null)
                this.kmers = new ProteinKmers(this.protein);
            // Compute the representation.
            retVal = repGenSet.findClosest(this.kmers);
        }
        return retVal;
    }

}
