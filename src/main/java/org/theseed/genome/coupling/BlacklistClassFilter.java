/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.coupling.FeatureClass.Result;
import org.theseed.io.TabbedLineReader;

/**
 * This filter removes classes found in a tab-delimited blacklist file.  The file is read in during the constructor
 * and applied to each genome.
 *
 * @author Bruce Parrello
 *
 */
public class BlacklistClassFilter extends ClassFilter {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(BlacklistClassFilter.class);
    /** set of blacklisted classes */
    private final Set<String> blackList;
    /** source blacklist file */
    private final File blackListFile;

    public BlacklistClassFilter(CouplesProcessor processor) throws IOException, ParseFailureException {
        super(processor);
        this.blackList = new HashSet<>();
        // Get the blacklist file.
        this.blackListFile = processor.getBlackListFile();
        if (blackListFile == null)
            throw new ParseFailureException("Blacklist file is required for filter type BLACKLIST.");
        if (! blackListFile.canRead())
            throw new FileNotFoundException("Blacklist file " + blackListFile + " is not found or unreadable.");
        else try (TabbedLineReader blackStream = new TabbedLineReader(blackListFile)) {
            for (TabbedLineReader.Line line : blackStream)
                this.blackList.add(line.get(0));
            log.info("{} blacklisted class IDs found in {}.", this.blackList.size(), blackListFile);
        }
    }

    @Override
    protected int apply(List<Result> results) {
        return this.removeClasses(this.blackList, results);
    }

    @Override
    protected String getName() {
        return String.format("BLACKLIST removing %d classes read from %s", this.blackList.size(), this.blackListFile);
    }

}
