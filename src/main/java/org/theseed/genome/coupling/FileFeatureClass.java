/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
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
public class FileFeatureClass extends FeatureClass {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FileFeatureClass.class);
    /** map of protein MD5s to family IDs */
    private Map<String, String> md5Map;
    /** map of family IDs to names */
    private Map<String, String> nameMap;

    /**
     * Initialize this classifier. We read in the family file and create the hashes.
     */
    public FileFeatureClass(FeatureClass.IParms processor) {
        // Create the maps.
        this.md5Map = new HashMap<String, String>();
        this.nameMap = new HashMap<String, String>();
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
                this.nameMap.put(famId, product);
            }
            log.info("{} records read. {} families for {} proteins.", inCount, this.nameMap.size(), this.md5Map.size());
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

    @Override
    public String getName(String classId) {
        return classId + "\t" + this.nameMap.getOrDefault(classId, "");
    }

    @Override
    public String getHeadings() {
        return("family_id1\tfamily_product1\tfamily_id2\tfamily_product2");
    }

    @Override
    public void cacheNames(Collection<String> classes) {
        // No need to cache, as all the names are permanently in memory.
    }

    @Override
    public Pair readPair(TabbedLineReader.Line line) {
        // The output line for protein families is of the form ID1 NAME1 ID2 NAME2.  We cache the names.
        String class1 = line.get(0);
        this.nameMap.put(class1, line.get(1));
        String class2 = line.get(2);
        this.nameMap.put(class2, line.get(3));
        return this.new Pair(class1, class2);
    }

}
