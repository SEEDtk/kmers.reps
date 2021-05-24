/**
 *
 */
package org.theseed.proteins;

/**
 * This enum is used to classify role counts when processing universal roles.  The count types are
 * missing, present exactly once, or over-present.
 *
 * @author Bruce Parrello
 *
 */
public enum RoleType {
    MISSING, GOOD, OVER_PRESENT;

    public static RoleType compute(int count) {
        RoleType retVal;
        switch (count) {
        case 0:
            retVal = MISSING;
            break;
        case 1:
            retVal = GOOD;
            break;
        default:
            retVal = OVER_PRESENT;
        }
        return retVal;
    }
}
