/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Set;

import org.theseed.genome.Genome;
import org.theseed.genome.TaxItem;


/**
 * This object produces a report on the taxonomy of each genome for a genome source.  The output is a simple
 * 3-column table:  genome ID, genome name, and taxonomy string.
 *
 * The constructor specifies the lowest rank to report; so, if the lowest rank is "genus", no species information
 * will be output, and if the lowest rank is "species", no strain information will be output.
 *
 * @author Bruce Parrello
 *
 */
public class RepTaxonReport {

    // FIELDS
    /** output writer */
    private PrintWriter writer;
    /** lowest rank to output */
    private Set<String> startRankSet;
    /** buffer for output line */
    private StringBuilder buffer;

    /**
     * Enumeration of supported minimum ranks.
     */
    public static enum MinRank {
        FAMILY {
            @Override
            public Set<String> stopSet() {
                return Set.of("family", "superfamily", "suborder", "order", "subclass", "class", "phylum");
            }
        }, GENUS {
            @Override
            public Set<String> stopSet() {
                return Set.of("genus", "subfamily", "family", "superfamily", "suborder", "order", "subclass", "class", "phylum");
            }
        }, SPECIES {
            @Override
            public Set<String> stopSet() {
                return Set.of("species", "genus", "subfamily", "family", "superfamily", "suborder", "order", "subclass", "class", "phylum");
            }
        };

        /**
         * @return the set of ranks which can possibly sit at the bottom of a lineage with the specified minimum rank
         */
        public abstract Set<String> stopSet();
    }

    /**
     * Construct a report on an output stream.
     *
     * @param outWriter		output print writer for the report
     * @param rank			lowest acceptable rank for the taxonomy string
     */
    public RepTaxonReport(PrintWriter outWriter, MinRank rank) {
        this.writer = outWriter;
        this.startRankSet = rank.stopSet();
        this.writer.println("genome_id\tgenome_name\ttaxonomy");
        this.buffer = new StringBuilder(500);
    }

    /**
     * Write out a data line for a genome.
     *
     * @param genome	genome to write
     */
    public void write(Genome genome) {
        // Create the ID and name columns.
        this.buffer.setLength(0);
        this.buffer.append(genome.getId());
        this.buffer.append('\t');
        this.buffer.append(genome.getName());
        this.buffer.append('\t');
        // Now loop through the lineage, building the taxonomy string.  First,
        // we find the minimum rank for this report.
        Iterator<TaxItem> taxIter = genome.taxonomy();
        TaxItem item = null;
        while (item == null && taxIter.hasNext()) {
            item = taxIter.next();
            if (! this.startRankSet.contains(item.getRank()))
                item = null;
        }
        // Now the item is either NULL or is the first taxonomic group we want to output.
        if (item == null)
            this.buffer.append(genome.getDomain());
        else {
            // Put in the smallest group.
            this.buffer.append(item.getName());
            while (taxIter.hasNext()) {
                this.buffer.append(", ");
                this.buffer.append(taxIter.next().getName());
            }
        }
        // Write the line.
        writer.println(this.buffer.toString());
    }
}
