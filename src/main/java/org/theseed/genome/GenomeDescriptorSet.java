/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.TreeSet;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.ParseFailureException;

/**
 * This class implements a set of genome descriptors.  Its primary function is to provide functionality for
 * finding the closest genome in the set.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeDescriptorSet extends TreeSet<GenomeDescriptor> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeDescriptorSet.class);
    /** similarity-based measures */
    public static int[] SIM_TYPES = new int[] { FinderType.SEED_SIMILARITY.ordinal(), FinderType.RNA_SIMILARITY.ordinal() };
    /** distance-based measures */
    public static int[] DIST_TYPES = new int[] { FinderType.SEED_DISTANCE.ordinal(), FinderType.RNA_DISTANCE.ordinal() };
    /** seed-based measures */
    public static int[] SEED_TYPES = new int[] { FinderType.SEED_SIMILARITY.ordinal(), FinderType.SEED_DISTANCE.ordinal() };
    /** serialization version ID */
    private static final long serialVersionUID = 8898552012420344740L;

    /**
     * Construct an empty genome descriptor set.
     */
    public GenomeDescriptorSet() {
        super();
    }

    /**
     * Construct a genome descriptor set from a four-column table.
     *
     * @param inFile	file containing the four-column table
     */
    public GenomeDescriptorSet(File inFile) throws IOException {
        try (TabbedLineReader reader = new TabbedLineReader(inFile)) {
            int idCol = reader.findField("genome_id");
            int nameCol = reader.findField("genome_name");
            int seedCol = reader.findField("seed_protein");
            int rnaCol = reader.findField("ssu_rna");
            for (TabbedLineReader.Line line : reader) {
                GenomeDescriptor newGenome = new GenomeDescriptor(line.get(idCol), line.get(nameCol),
                        line.get(seedCol), line.get(rnaCol));
                this.add(newGenome);
            }
        }
        log.info("{} genomes read from {}.", this.size(), inFile);
    }

    /**
     * Add a genome to this set.
     *
     * @param genome	genome to add
     *
     * @throws ParseFailureException
     */
    public void add(Genome genome) throws ParseFailureException {
        GenomeDescriptor newGenome = new GenomeDescriptor(genome);
        this.add(newGenome);
    }

    /**
     * This enum determines the types of comparisons that can be made between genomes.  All of these are converted into
     * a "proximity" rating.  A higher proximity means the genomes are closer.
     *
     * @author Bruce Parrello
     *
     */
    public static enum FinderType {
        SEED_SIMILARITY {
            @Override
            public double getProximity(GenomeDescriptor testGenome, GenomeDescriptor refGenome) {
                return (double) testGenome.getSeedSim(refGenome);
            }

            @Override
            public String label() {
                return "seed_sim";
            }
        }, RNA_SIMILARITY {
            @Override
            public double getProximity(GenomeDescriptor testGenome, GenomeDescriptor refGenome) {
                return (double) testGenome.getRnaSim(refGenome);
            }

            @Override
            public String label() {
                return "rna_sim";
            }
        }, SEED_DISTANCE {
            @Override
            public double getProximity(GenomeDescriptor testGenome, GenomeDescriptor refGenome) {
                return 1.0 - testGenome.getSeedDistance(refGenome);
            }

            @Override
            public String label() {
                return "seed_distance";
            }
        }, RNA_DISTANCE {
            @Override
            public double getProximity(GenomeDescriptor testGenome, GenomeDescriptor refGenome) {
                return 1.0 - testGenome.getRnaDistance(refGenome);
            }

            @Override
            public String label() {
                return "rna_distance";
            }
        };

        /**
         * @return the proximity rating of a test genome and a reference genome
         *
         * @param testGenome	test genome's descriptor
         * @param refGenome		reference genome's descriptor
         */
        public abstract double getProximity(GenomeDescriptor testGenome, GenomeDescriptor refGenome);

        /**
         * @return a descriptor string for this type, suitable for column headers
         */
        public abstract String label();
    }

    /**
     * This class represents the results of an attempt to find the closest genome in the set.
     */
    public static class CloseGenome implements Comparable<CloseGenome> {

        private GenomeDescriptor genome;
        private double proximity;

        /**
         * Construct an empty close-genome result.
         */
        protected CloseGenome() {
            this.genome = null;
            this.proximity = 0.0;
        }

        /**
         * The closest genome (highest proximity) always sorts first.
         */
        @Override
        public int compareTo(CloseGenome o) {
            int retVal = Double.compare(o.proximity, this.proximity);
            if (retVal == 0)
                retVal = this.genome.compareTo(o.genome);
            return retVal;
        }

        /**
         * Update this result with the new genome if it's closer.
         *
         * @param genome	descriptor of new genome
         * @param type		type of proximity measure
         *
         * @return TRUE if the new genome is closer
         */
        public boolean check(GenomeDescriptor testGenome, GenomeDescriptor newRefGenome, FinderType type) {
            double newProx = type.getProximity(testGenome, newRefGenome);
            boolean retVal = false;
            if (newProx > this.proximity) {
                retVal = true;
                this.proximity = newProx;
                this.genome = newRefGenome;
            }
            return retVal;
        }

        /**
         * @return a tab-delimited output string for this result
         */
        public String output() {
            StringBuilder retVal = new StringBuilder(40);
            if (this.genome != null) {
                retVal.append(this.genome.getId());
            }
            retVal.append("\t");
            retVal.append(String.format("%8.4f", this.proximity));
            return retVal.toString();
        }

        /**
         * @return TRUE if this result has the same genome as another result
         *
         * @param other		other result to compare
         */
        public boolean isSameGenome(CloseGenome other) {
            boolean retVal;
            if (this.genome == null)
                retVal = (other.genome == null);
            else if (other.genome == null)
                retVal = false;
            else
                retVal = (this.genome.getId().contentEquals(other.genome.getId()));
            return retVal;
        }

        /**
         * @return TRUE if all the close genomes in a collection are the same, else FALSE
         *
         * @param results		a list of the results for all the types
         * @param positions		an array of the relevant positions in the list
         */
        public static boolean test(List<CloseGenome> results, int[] positions) {
            boolean retVal = true;
            CloseGenome result0 = results.get(positions[0]);
            for (int i = 1; i < positions.length && retVal; i++) {
                retVal = result0.isSameGenome(results.get(positions[i]));
            }
            return retVal;
        }

        /**
         * @return the proximity
         */
        public double getProximity() {
            return this.proximity;
        }

        /**
         * @return the genome ID of the close genome
         */
        public String getGenomeId() {
            String retVal;
            if (this.genome == null)
                retVal = "";
            else
                retVal = this.genome.getId();
            return retVal;
        }

        /**
         * @return the name of the close genome
         */
        public String getGenomeName() {
            String retVal;
            if (this.genome == null)
                retVal = "(none)";
            else
                retVal = this.genome.getName();
            return retVal;
        }

        /**
         * @return the descriptor of the close genome
         */
        public GenomeDescriptor getGenome() {
            return this.genome;
        }

    }

    /**
     * Find the closest genome to the specified testing genome according to the specified finder type.
     *
     * @param genome	descriptor for the testing genome
     * @param type		type of proximity criterion to use
     *
     * @return a close-genome object describing the result
     */
    public CloseGenome findClosest(GenomeDescriptor testGenome, FinderType type) {
        // Start with no result.
        CloseGenome retVal = new CloseGenome();
        // Loop through the set, finding the closest.
        for (GenomeDescriptor refGenome : this)
            retVal.check(testGenome, refGenome, type);
        return retVal;
    }

}
