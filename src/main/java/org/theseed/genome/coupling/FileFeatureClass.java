/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.HashMap;
import java.util.Map;

import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;

/**
 * This classification system reads in a file of MD5s. Each is associated with an ID and a name. These IDs and
 * names are used as the protein family IDs and names assigned to features. Note that since each MD5 can only
 * occur once, only one class per feature is returned.
 *
 * @author Bruce Parrello
 *
 */
public class FileFeatureClass extends MapFeatureClass {

    // FIELDS
    /** map of protein MD5s to family IDs */
    private Map<String, String> md5Map;
    /**
     * Initialize this classifier. We read in the family file and create the hashes.
     */
    public FileFeatureClass(FeatureClass.IParms processor) {
        super();
        // Create the maps.
        this.md5Map = new HashMap<String, String>();
        File famFile = processor.getFamFile();
        log.info("Reading family data from {}.", famFile);
        try (TabbedLineReader inStream = new TabbedLineReader(famFile)) {
            // Find the input columns we need.
            int idIdx = inStream.findField("fam_id");
            int nameIdx = inStream.findField("product");
            int md5Idx = inStream.findField("md5");
            // Loop through the file, filling the maps.
            int inCount = 0;
            for (var line : inStream) {
                inCount++;
                String famId = line.get(idIdx);
                String product = line.get(nameIdx);
                String md5 = line.get(md5Idx);
                this.md5Map.put(md5, famId);
                this.put(famId, product);
            }
            log.info("{} records read. {} families for {} proteins.", inCount, this.size(), this.md5Map.size());
        } catch (IOException e) {
            // Uncheck the exception to conform to the API.
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public Result getClasses(Feature feat) {
        Result retVal = new Result(feat);
        // Find the protein MD5 and compute the family.
        String md5 = feat.getMD5();
        String family = this.md5Map.get(md5);
        if (family != null)
            retVal.add(family);
        return retVal;
    }

}
