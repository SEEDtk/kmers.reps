/**
 *
 */
package org.theseed.genome.coupling;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.theseed.genome.Feature;
import org.theseed.p3api.Connection;
import org.theseed.p3api.Connection.Table;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This is a classification system based on protein families.  For protein families, the name is loaded from PATRIC.
 *
 * @author Bruce Parrello
 *
 */
public class FamilyFeatureClass extends FeatureClass {

    // FIELDS
    /** connection to PATRIC */
    private Connection p3;
    /** map of family IDs to names */
    private Map<String, String> nameMap;

    /**
     * Construct the family feature class.
     */
    public FamilyFeatureClass() {
        this.p3 = new Connection();
        this.nameMap = new HashMap<String, String>(5000);
    }

    @Override
    public Result getClasses(Feature feat) {
        Result retVal = null;
        String family = feat.getPgfam();
        if (family != null) {
            retVal = new Result(feat);
            retVal.add(family);
        }
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
                this.nameMap.put(familyDatum.getKey(), Connection.getString(familyDatum.getValue(), "family_product"));
        }
    }

}