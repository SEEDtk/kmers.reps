/**
 *
 */
package org.theseed.proteins;

/**
 * This class converts a protein into an array of booleans indicating which 3mers are present or
 * absent in the protein.  The resulting array can be converted into a feature for input into a
 * learning system.
 *
 * @author Bruce Parrello
 *
 */
public final class ProteinEncoder {

    /** array for converting amino acid codes to digits */
    private static final String AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY";

    /** constant representing a null array index */
    private static final int ANULL = -1;

    /**
     * @return the header string for a file containing encoded proteins
     */
    public static String header() {
        // The result will be built in here.
        StringBuilder retVal = new StringBuilder(AMINO_ACIDS.length() * 4);
        for (int i = 0; i < AMINO_ACIDS.length(); i++) {
            for (int j = 0; j < AMINO_ACIDS.length(); j++) {
                for (int k = 0; k < AMINO_ACIDS.length(); k++) {
                    retVal.append(AMINO_ACIDS.charAt(i));
                    retVal.append(AMINO_ACIDS.charAt(j));
                    retVal.append(AMINO_ACIDS.charAt(k));
                    retVal.append('\t');
                }
            }
        }
        // Return everything but the extra tab at the end.
        return retVal.substring(0, AMINO_ACIDS.length() * 4 - 1);
    }

    /**
     * @return the index number corresponding to a 3-mer, or ANULL if the 3-mer is invalid
     *
     * @param prot		the protein sequence of interest
     * @param pos		the position (1-based) of the 3-mer in the protein
     */
    public static int getIdx(String prot, int pos) {
        int retVal = 0;
        int last = pos + 1;
        for (int i = pos - 1; i <= last && retVal >= 0; i++) {
            int idx = AMINO_ACIDS.indexOf(prot.charAt(i));
            if (idx < 0) {
                retVal = ANULL;
            } else {
                retVal = retVal * AMINO_ACIDS.length() + idx;
            }
        }
        return retVal;
    }

}
