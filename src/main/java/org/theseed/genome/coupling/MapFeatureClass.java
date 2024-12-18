/**
 *
 */
package org.theseed.genome.coupling;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.theseed.io.TabbedLineReader;

/**
 * The map-feature class is a common superclass for family schemes where the family names are read from a file.
 * @author Bruce Parrello
 *
 */
public abstract class MapFeatureClass extends FeatureClass {

    /** map of family IDs to names */
    private Map<String, String> nameMap;

    /**
     * Construct a map-feature class.
     */
    public MapFeatureClass() {
        super();
        // Create an empty name map.
        this.nameMap = new HashMap<String, String>();
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

    /**
     * This is a service method for putting names in the map.
     *
     * @param id		family ID
     * @param name		family name
     */
    protected void put(String id, String name) {
        this.nameMap.put(id, name);
    }

    /**
     * @return the name for a family
     *
     * @param id		family ID
     */
    protected String get(String id) {
        return this.nameMap.getOrDefault(id, "");
    }

    /**
     * @return the number of families in the name map
     */
    public int size() {
        return this.nameMap.size();
    }

}
