/**
 *
 */
package org.theseed.genome;

import org.theseed.proteins.kmers.reps.SeqTableProcessor;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.ParseFailureException;

/**
 * This class contains the genome ID, name, seed protein, and SSU rRNA sequence for a genome.  This is the basic
 * unit of information for a four-column table output by the seqTable command.  (It should be noted that the actual
 * table contains 5 columns.)
 *
 * @author Bruce Parrello
 *
 */
public class GenomeDescriptor implements Comparable<GenomeDescriptor> {

    // FIELDS
    /** genome ID */
    private String genomeId;
    /** genome name */
    private String name;
    /** seed protein kmers */
    private ProteinKmers seedKmers;
    /** SSU rRNA kmers */
    private DnaKmers rnaKmers;

    /**
     * Construct a descriptor from a genome.
     *
     * @param genome	genome from which to build the descriptor
     *
     * @throws ParseFailureException	if a sequence is missing
     */
    public GenomeDescriptor(Genome genome) throws ParseFailureException {
        this.genomeId = genome.getId();
        this.name = genome.getName();
        // Find the seed protein and build kmers from its protein translation.
        String protId = SeqTableProcessor.findSeed(genome);
        if (protId == null)
            throw new ParseFailureException("Genome " + this.genomeId + "(" + this.name + ") is missing a seed protein.");
        Feature seed = genome.getFeature(protId);
        this.seedKmers = new ProteinKmers(seed.getProteinTranslation());
        // Find the SSU rRNA and build kmers from its DNA.
        String rna = genome.getSsuRRna();
        if (rna.isEmpty())
            throw new ParseFailureException("Genome " + this.genomeId + "(" + this.name + ") is missing an SSU rRNA sequence.");
        this.rnaKmers = new DnaKmers(genome.getSsuRRna());
    }

    /**
     * Construct a descriptor from the basic information.
     *
     * @param id		ID of the genome
     * @param name		name of the genome
     * @param seedProt	seed protein string
     * @param ssuRna	SSU rRNA DNA sequence
     */
    public GenomeDescriptor(String id, String name, String seedProt, String ssuRna) {
        this.genomeId = id;
        this.name = name;
        this.seedKmers = new ProteinKmers(seedProt);
        this.rnaKmers = new DnaKmers(ssuRna);
    }

    /**
     * Compute the seed protein similarity to another genome.
     *
     * @param other		other descriptor to check
     *
     * @return the number of seed protein kmers in common
     */
    public int getSeedSim(GenomeDescriptor other) {
        return this.seedKmers.similarity(other.seedKmers);
    }

    /**
     * Compute the SSU rRNA similarity to another genome.
     *
     * @param other		other descriptor to check
     *
     * @return the number of SSU rRNA kmers in common
     */
    public int getRnaSim(GenomeDescriptor other) {
        return this.rnaKmers.similarity(other.rnaKmers);
    }

    /**
     * Compute the seed protein distance to another genome.
     *
     * @param other		other descriptor to check
     *
     * @return the kmer distance between the seed proteins of the genomes
     */
    public double getSeedDistance(GenomeDescriptor other) {
        return this.seedKmers.distance(other.seedKmers);
    }

    /**
     * Compute the SSU rRNA distance to another genome.
     *
     * @param other		other descriptor to check
     *
     * @return the kmer distance between the SSU rRNA sequences of the genomes
     */
    public double getRnaDistance(GenomeDescriptor other) {
        return this.rnaKmers.distance(other.rnaKmers);
    }

    /**
     * @return the genome ID
     */
    public String getId() {
        return this.genomeId;
    }

    /**
     * @return the genome name
     */
    public String getName() {
        return this.name;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.genomeId == null) ? 0 : this.genomeId.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof GenomeDescriptor)) {
            return false;
        }
        GenomeDescriptor other = (GenomeDescriptor) obj;
        if (this.genomeId == null) {
            if (other.genomeId != null) {
                return false;
            }
        } else if (!this.genomeId.equals(other.genomeId)) {
            return false;
        }
        return true;
    }

    /**
     * We sort by name and then ID.  Ideally, we want to sort by taxonomy, but we don't have a way to
     * get that cheaply from a standard four-column table file.
     */
    @Override
    public int compareTo(GenomeDescriptor o) {
        int retVal = this.name.compareTo(o.name);
        if (retVal == 0)
            retVal = this.genomeId.compareTo(o.genomeId);
        return retVal;
    }

}
