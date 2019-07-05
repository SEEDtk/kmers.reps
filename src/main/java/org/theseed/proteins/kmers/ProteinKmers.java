/**
 *
 */
package org.theseed.proteins.kmers;

import java.util.HashSet;

/**
 * This package manages the set of protein kmers in a protein string.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinKmers {

    /** current kmer size */
    private static int K = 10;

    // FIELDS
    /** initial protein string */
    private String protein;
    /** set of kmers in the protein */
    private HashSet<String> kmerSet;

    /**
     * Generate a protein kmer set for a specified protein.
     */
    public ProteinKmers(String protein) {
        this.protein = protein;
        int n = protein.length() - K;
        this.kmerSet = new HashSet<String>(n);
        for (int i = 0; i <= n; i++) {
            kmerSet.add(protein.substring(i, i + K));
        }
    }

    /**
     * Specify a new global protein kmer size.
     *
     * @param kSize	proposed new kmer size
     */
    public static void setKmerSize(int kSize) {
        K = kSize;
    }

    /**
     * @return the current protein kmer size
     */
    public static int kmerSize() {
        return K;
    }

    /**
     * @return the number of kmers in common between two proteins
     *
     * @param other		the protein-kmers object for the other protein
     */
    public int similarity(ProteinKmers other) {
        int retVal = 0;
        for (String kmer : other.kmerSet) {
            if (this.kmerSet.contains(kmer)) {
                retVal++;
            }
        }
        return retVal;
    }

    /**
     * @return the Jaccard distance between two proteins
     *
     * @param other		the protein-kmers object for the other protein
     */
    public double distance(ProteinKmers other) {
        double retVal = 0;
        double similarity = this.similarity(other);
        if (similarity > 0) {
            double union = (this.kmerSet.size() + other.kmerSet.size()) - similarity;
            retVal = 1.0 - similarity / union;
        }
        return retVal;
    }

    /**
     * @return the protein
     */
    public String getProtein() {
        return this.protein;
    }

    /**
     * @return the kmer count
     */
    public int size() {
        return this.kmerSet.size();
    }



}
