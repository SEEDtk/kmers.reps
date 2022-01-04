/**
 *
 */
package org.theseed.sequence;

import java.util.Comparator;

import org.theseed.genome.Genome;
import org.theseed.reports.NaturalSort;

/**
 * This object describes a genome representative that uses an SSU rather than a seed protein
 * sequence.  At some point, if we get serious about SSU analysis, we will want to unify the two
 * structures.
 *
 * The object is identified by genome ID, so the hash code and the comparison methods only look at
 * the genome ID.
 *
 * @author Bruce Parrello
 *
 */
public class SsuGenome extends DnaKmers implements Comparable<SsuGenome> {

    // FIELDS
    /** genome ID */
    private String genomeId;
    /** genome name */
    private String name;
    /** natural sorter for genome IDs */
    private static final Comparator<String> NATURAL_SORT = new NaturalSort();

    /**
     * Create an SSU representation object for a specified genome name and ID.
     *
     * @param genome_id		ID of the genome
     * @param gName			name of the genome
     * @param dna			SSU sequence to use
     */
    public SsuGenome(String genome_id, String gName, String dna) {
        super(dna);
        this.genomeId = genome_id;
        this.name = gName;
    }

    /**
     * Create an SSU object for a particular genome.
     *
     * @param genome	genome of interest
     *
     * @return the SSU object created, or NULL if the genome has no SSU
     */
    public static SsuGenome create(Genome genome) {
        SsuGenome retVal = null;
        String ssu = genome.getSsuRRna();
        if (! ssu.isEmpty()) {
            // Here we have a valid SSU, so we can create the object.
            retVal = new SsuGenome(genome.getId(), genome.getName(), ssu);
        }
        return retVal;
    }

    /**
     * @return the ID of the associated genome
     */
    public String getGenomeId() {
        return this.genomeId;
    }

    /**
     * @return the name of the associated genome
     */
    public String getName() {
        return this.name;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = super.hashCode();
        result = prime * result + ((this.genomeId == null) ? 0 : this.genomeId.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!super.equals(obj)) {
            return false;
        }
        if (!(obj instanceof SsuGenome)) {
            return false;
        }
        SsuGenome other = (SsuGenome) obj;
        if (this.genomeId == null) {
            if (other.genomeId != null) {
                return false;
            }
        } else if (!this.genomeId.equals(other.genomeId)) {
            return false;
        }
        return true;
    }

    @Override
    public int compareTo(SsuGenome o) {
        return NATURAL_SORT.compare(this.genomeId, o.genomeId);
    }

}
