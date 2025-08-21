/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;

/**
 * This classification method uses local protein families developed during protein family testing. The family
 * names are taken from a definition file that has the family ID in column 1 and the name in column 2.
 *
 * @author Bruce Parrello
 *
 */
public class TestLocalFamilyFeatureClass extends MapFeatureClass {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(TestLocalFamilyFeatureClass.class);

    public TestLocalFamilyFeatureClass(FeatureClass.IParms processor) {
        super();
        // Fill the name map from the input file.
        File famFile = processor.getFamFile();
        log.info("Reading family data from {}.", famFile);
        try (TabbedLineReader inStream = new TabbedLineReader(famFile)) {
            for (var line : inStream) {
                String famId = line.get(0);
                String name = line.get(1);
                this.put(famId, name);
            }
            log.info("{} family names read from {}.", this.size(), famFile);
        } catch (IOException e) {
            // Uncheck the exception to conform to the API.
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public Result getClasses(Feature feat) {
        Result retVal = new Result(feat);
        String family = feat.getPlfam();
        if (! StringUtils.isBlank(family))
            retVal.add(family);
        return retVal;
    }

}
