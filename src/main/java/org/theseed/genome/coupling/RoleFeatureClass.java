/**
 *
 */
package org.theseed.genome.coupling;

import java.util.Collection;

import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader.Line;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * This feature classifier returns all the roles in the feature.
 *
 * @author Bruce Parrello
 *
 */
public class RoleFeatureClass extends FeatureClass {

    // FIELDS
    /** map of roles */
    private RoleMap roles;

    /**
     * Construct a feature role classifier.
     */
    public RoleFeatureClass() {
        // We build an empty role map and add non-hypotheticals as we go.
        this.roles = new RoleMap();
    }

    @Override
    public Result getClasses(Feature feat) {
        // Create a prototype result.
        Result retVal = new Result(feat);
        // Get the feature's function and separate out the roles.
        String function = feat.getFunction();
        if (function != null  && ! function.isEmpty()) {
            String[] roles = Feature.rolesOfFunction(function);
            // Convert the role names to IDs.  If the role is new, we add it to the map.
            for (String role : roles) {
                if (! Feature.isHypothetical(role)) {
                    String roleId = this.roles.findOrInsert(role).getId();
                    retVal.add(roleId);
                }
            }
        }
        // Return NULL if no roles were found.
        if (! retVal.isGood()) retVal = null;
        return retVal;
    }

    @Override
    public String getName(String classId) {
        String retVal = this.roles.getName(classId);
        if (retVal == null) retVal = "";
        return retVal;
    }

    @Override
    public String getHeadings() {
        return("role1\trole2");
    }

    @Override
    public void cacheNames(Collection<String> classes) {
        // No name cache is needed for roles, since they exist in the map.
    }

    @Override
    public Pair readPair(Line line) {
        // The output line consists of NAME1 NAME2.  We need the IDs.
        Role role1 = this.roles.findOrInsert(line.get(0));
        Role role2 = this.roles.findOrInsert(line.get(1));
        return this.new Pair(role1.getId(), role2.getId());
    }

}
