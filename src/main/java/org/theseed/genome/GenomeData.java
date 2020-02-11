/**
 *
 */
package org.theseed.genome;

/**
 * This is a very simple class that represents the data we need about a genome-- genus, species, and seed protein.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeData {

    public static final int GENUS = 0;
    public static final int SPECIES = 1;
    // FIELDS
    private String taxon[];
    private String seedProtein;

    /**
     * Create a new genome data element.
     */
    public GenomeData() {
        this.taxon = new String[] { null, null };
        this.seedProtein = null;
    }

    /**
     * @return the seed protein sequence
     */
    public String getSeedProtein() {
        return seedProtein;
    }

    /**
     * Store the seed protein.
     *
     * @param seedProtein the seed protein string
     */
    public GenomeData setSeedProtein(String seedProtein) {
        this.seedProtein = seedProtein;
        return this;
    }

    /**
     * Store the genus and species.
     *
     * @param genus the genus to set
     */
    public GenomeData setTaxonomy(String genus, String species) {
        this.taxon[GENUS] = genus;
        this.taxon[SPECIES] = species;
        return this;
    }

    /**
     * @return the specified taxon ID
     *
     * @param idx	type of taxon ID
     */
    public String getTaxon(int idx) {
        return taxon[idx];
    }

}
