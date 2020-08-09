/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.HashSet;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.coupling.FeatureClass.Pair;
import org.theseed.genome.coupling.FeatureClass.PairData;
import org.theseed.io.TabbedLineReader;


/**
 * This filter accepts as input a whitelist of class IDs.  Only pairs for which at least one class is included
 * in the whitelist will be considered significant.
 *
 * @author Bruce Parrello
 *
 */
public class WhitePairFilter extends PairFilter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(WhitePairFilter.class);
    /** set of acceptable classes */
    private Set<String> whiteSet;

    public WhitePairFilter(CouplesProcessor processor) {
        super(processor);
        // Get the input file.
        File whiteFile = processor.getWhiteGroupFile();
        if (whiteFile == null)
            throw new IllegalArgumentException("Whitelist file required for pair-filtering of type WHITELIST.");
        try (TabbedLineReader whiteStream = new TabbedLineReader(whiteFile)) {
            this.whiteSet = new HashSet<String>();
            for (TabbedLineReader.Line line : whiteStream)
                this.whiteSet.add(line.get(0));
            log.info("{} class identifiers in whitelist.", this.whiteSet.size());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public boolean isSignificant(Pair pair, PairData pairData) {
        return (this.whiteSet.contains(pair.getClass1()) || this.whiteSet.contains(pair.getClass2()));
    }

    @Override
    protected String getName() {
        return "WHITELIST";
    }

}
