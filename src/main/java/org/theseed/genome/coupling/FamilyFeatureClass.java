/**
 *
 */
package org.theseed.genome.coupling;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader.Line;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Connection.Table;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This is a classification system based on protein families.  For protein families, the name is loaded from PATRIC, or
 * optionally cached from a file.
 *
 * @author Bruce Parrello
 *
 */
public class FamilyFeatureClass extends FeatureClass {

    // FIELDS
    /** logging facility */
    protected Logger log = LoggerFactory.getLogger(FamilyFeatureClass.class);
    /** connection to PATRIC */
    private P3Connection p3;
    /** map of family IDs to names */
    private Map<String, String> nameMap;

    /**
     * Construct the family feature class.
     *
     * @param processor		controlling command processor
     */
    public FamilyFeatureClass(IParms processor) {
        this.p3 = new P3Connection();
        this.nameMap = new HashMap<String, String>(5000);
    }

    @Override
    public Result getClasses(Feature feat) {
        Result retVal = new Result(feat);
        String family = feat.getPgfam();
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
        // Gather the unknown classes in the incoming list.
        Set<String> batch = new HashSet<String>(classes.size());
        for (String classId : classes) {
            if (! this.nameMap.containsKey(classId))
                batch.add(classId);
        }
        if (batch.size() > 0) {
            Map<String, JsonObject> familyData = this.p3.getRecords(Table.FAMILY, batch, "family_product");
            for (Map.Entry<String, JsonObject> familyDatum : familyData.entrySet())
                this.nameMap.put(familyDatum.getKey(), KeyBuffer.getString(familyDatum.getValue(), "family_product"));
        }
    }

    @Override
    public Pair readPair(Line line) {
        // The output line for protein families is of the form ID1 NAME1 ID2 NAME2.  We cache the names.
        String class1 = line.get(0);
        this.nameMap.put(class1, line.get(1));
        String class2 = line.get(2);
        this.nameMap.put(class2, line.get(3));
        return this.new Pair(class1, class2);
    }

    /**
     * Read IDs and names from a FASTA file and cache them in this object.  This method allows
     * specification of families not in PATRIC.
     *
     * @param fastaFile		FASTA file to read
     */
    public void cacheNames(File fastaFile) {
        log.info("Reading family names from {}.", fastaFile);
        try (FastaInputStream fastaStream = new FastaInputStream(fastaFile)) {
            int count = 0;
            for (Sequence seq : fastaStream) {
                this.nameMap.put(seq.getLabel(), seq.getComment());
                count++;
            }
            log.info("{} family names read from file.", count);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

}
