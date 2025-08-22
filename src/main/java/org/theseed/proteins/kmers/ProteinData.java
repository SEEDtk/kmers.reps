/**
 *
 */
package org.theseed.proteins.kmers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.text.TextStringBuilder;
import org.theseed.genome.Genome;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.ProteinKmers;

/**
 * This class describes the seed protein for a genome.  It contains the DNA and amino acid sequences,
 * various information about the genome, and the SSU rRNA sequence.  It is sorted by score within rating.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinData implements Comparable<ProteinData> {

    /** feature ID */
    private String fid;
    /** genome ID */
    private final String genomeId;
    /** genome name */
    private final String genomeName;
    /** DNA sequence */
    private String dna;
    /** amino acid sequence */
    private String protein;
    /** taxonomic genus */
    private String genus;
    /** taxonomic species */
    private String species;
    /** quality score */
    private final double score;
    /** rating */
    private Rating rating;
    /** domain name */
    private String domain;
    /** genetic code */
    private final int geneticCode;
    /** proteinkmers object for computing representation */
    private ProteinKmers kmers;
    /** sequence for the best SSU rRNA */
    private String ssuSequence;
    /** representative genome information for each score level */
    private final Map<Integer, RepGenomeDb.Representation> representation;

    /**
     * This enum gives the rating of a protein data object.  The ratings are ordered from
     * best to worst.
     */
    public static enum Rating {
        /** flagged as an NCBI reference genome */
        NCBI_REF,
        /** flagged as an NCBI representative genome */
        NCBI_REP,
        /** good PATRIC genome */
        NORMAL,
        /** short SSU */
        SHORT_SSU,
        /** single SSU */
        SINGLE_SSU,
        /** bad or questionable SSU */
        BAD_SSU;

    }

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
        this.rating = Rating.NORMAL;
        this.representation = new HashMap<>();
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
     * @return the protein sequence
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

    /**
     * Find the closest representative of this genome in a repgen set.
     *
     * @param repGenSet			repgen set to check
     *
     * @return the representation information for this genome.
     */
    public RepGenomeDb.Representation getRepresentation(RepGenomeDb repGenSet) {
        int sim = repGenSet.getThreshold();
        RepGenomeDb.Representation retVal = this.representation.get(sim);
        // If we don't already have a representative, compute one.
        if (retVal == null) {
            // Insure we have the protein kmers.
            if (this.kmers == null)
                this.kmers = new ProteinKmers(this.protein);
            if (retVal == null) {
                // We need a new representative. Compute it here.
                retVal = repGenSet.findClosest(this.kmers);
            }
        }
        return retVal;
    }

    /**
     * Update the taxonomic information for this genome.
     *
     * @param domain	domain of the genome
     * @param genus		genus of the genome
     * @param species	specieis of the genome
     */
    public void setTaxonomy(String domain, String genus, String species) {
        this.domain = domain;
        this.genus = genus;
        this.species = species;
    }

    /**
     * Store the representation information for this genome at the specified similarity threshold level.
     *
     * @param repDb			relevant representative-genome database
     * @param repId			ID of the representing genome
     * @param score			similarity score between this genome and the representative
     * @param distance		distance between this genome and the representative
     */
    public void setRepresentation(RepGenomeDb repDb, String repId, int score, double distance) {
        RepGenomeDb.Representation newRep = repDb.new Representation(repId, score, distance);
        this.setRepresentation(repDb, newRep);
    }

    /**
     * Store the representation information for this genome at the specified similarity threshold level.
     *
     * @param repDb			relevant representative-genome database
     * @param rep			representation to store
     */
    public void setRepresentation(RepGenomeDb repDb, RepGenomeDb.Representation rep) {
        this.representation.put(repDb.getThreshold(), rep);
    }

    @Override
    public int compareTo(ProteinData other) {
        // The best rating sorts first, then the highest score, and finally the genome ID.
        int retVal = this.rating.compareTo(other.rating);
        if (retVal == 0) {
            retVal = Double.compare(other.score, this.score);
            if (retVal == 0)
                retVal = this.genomeId.compareTo(other.genomeId);
        }
        return retVal;
    }

    /**
     * @return TRUE if this genome's SSU rRNA is questionable
     */
    public boolean isQuestionable() {
        return this.rating == Rating.BAD_SSU;
    }
    /**
     * @return the SSU rRNA sequence for this genome
     */
    public String getSsuSequence() {
        return this.ssuSequence;
    }

    /**
     * Specify a new genome rating.
     * @param rating 	the new genome rating to store
     */
    public void setRating(Rating rating) {
        this.rating = rating;
    }

    /**
     * Update the SSU rRNA sequence.  We keep the longest, and we mark the genome as bad if it
     * has fewer than two SSUs of sufficient length or the SSUs are too far apart.
     *
     * @param rnas		list of SSU sequences for this genome
     *
     * @return the final rating
     */
    public Rating setSsuSequence(List<String> rnas) {
        // Start with an empty SSU sequence.
        this.ssuSequence = "";
        // Create a place to stash the good sequences.
        List<DnaKmers> goodSeqs = new ArrayList<>(rnas.size());
        // Count the sequences of acceptable length.
        int longSeqs = 0;
        // Check all the sequences.
        for (String rna : rnas) {
            // Only proceeed if the sequence has no long ambiguity run.
            if (Genome.isValidSsuRRna(rna)) {
                // Insure we keep the longest.
                if (rna.length() > this.ssuSequence.length())
                    this.ssuSequence = rna;
                // Save the good ones.
                if (rna.length() >= ProteinDataFactory.USEFUL_SSU_LEN) {
                    goodSeqs.add(new DnaKmers(rna));
                    if (rna.length() >= ProteinDataFactory.MIN_SSU_LEN)
                        longSeqs++;
                }
            }
        }
        if (goodSeqs.size() < 1) {
            // Here all the SSUs are fragments.  This counts as a missing SSU.
            this.rating = Rating.BAD_SSU;
        } else if (goodSeqs.size() < 2) {
            // Here there are insufficient long SSUs to validate them.  If this is an NCBI
            // genome and the single SSU is long, we trust it automatically.  If it is not
            // NCBI, we mark it SINGLE for possible reclamation later.  If it is NCBI but its SSU
            // is short, we mark it short, which means it automatically gets lower priority
            // than genomes with a long SSU that are NCBI or independently verified.
            if (this.rating.compareTo(Rating.NORMAL) >= 0)
                this.rating = Rating.SINGLE_SSU;
            else if (longSeqs < 1)
                this.rating = Rating.SHORT_SSU;
        } else {
            // Compare all the good sequences.
            boolean okFlag = true;
            for (int i = 0; okFlag && i < goodSeqs.size(); i++) {
                DnaKmers ssuKmers = goodSeqs.get(i);
                for (int j = i + 1; okFlag && j < goodSeqs.size(); j++) {
                    double dist = ssuKmers.distance(goodSeqs.get(j));
                    okFlag = (dist <= ProteinDataFactory.MAX_SSU_DISTANCE);
                }
            }
            if (! okFlag)
                this.rating = Rating.BAD_SSU;
            else if (longSeqs < 1)
                this.rating = Rating.SHORT_SSU;
        }
        return this.rating;
    }

    /**
     * @return the rating of this genome
     */
    public Rating getRating() {
        return this.rating;
    }

    /**
     * @return the representative genome for the specified level
     *
     * @param level		representation threshold
     */
    public String getRepGenome(int level) {
        RepGenomeDb.Representation rep = this.representation.get(level);
        String retVal;
        if (rep == null)
            retVal = "";
        else
            retVal = rep.getGenomeId();
        return retVal;
    }

    /**
     * @return the header for a ProteinData report
     *
     * @param parent	parent protein data factory being written
     */
    public static String getHeader(ProteinDataFactory parent) {
        String retVal = "genome_id\tgenome_name\tscore\trating\tgc\tdomain\tgenus\tspecies\t"
                + Arrays.stream(parent.getRepLevels()).mapToObj(i -> Integer.toString(i)).collect(Collectors.joining("\t"));
        return retVal;
    }

    /**
     * @return a report data line for this protein
     *
     * @param parent	parent protein data factory being written
     */
    public String getLine(ProteinDataFactory parent) {
        TextStringBuilder retVal = new TextStringBuilder(100);
        // Store the fixed fields.
        retVal.append("%s\t%s\t%6.2f\t%s\t%s\t%s\t%s", this.genomeId, this.genomeName, this.score, this.rating.toString(),
                this.domain, this.genus, this.species);
        // Add the representative genomes.
        Arrays.stream(parent.getRepLevels()).mapToObj(i -> this.getRepGenome(i)).forEach(x -> retVal.append('\t').append(x));
        return retVal.toString();
    }

}
