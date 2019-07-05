/**
 *
 */
package org.theseed.proteins;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.nio.charset.StandardCharsets.*;

import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.codec.binary.Base64;
/**
 * This class represents a role.  It stores the role ID, the role description, and
 * its checksum.  A static method is provided to parse the checksum from the
 * role description.
 *
 * @author Bruce Parrello
 *
 */
public class Role implements Comparable<Role> {

    // FIELDS
    String roleId;
    String checksum;
    String description;

    // ROLE-PARSING PATTERNS
    static private final Pattern EC_PATTERN = Pattern.compile("(.+?)\\s*\\(\\s*E\\.?C\\.?(?:\\s+|:)(\\d\\.(?:\\d+|-)\\.(?:\\d+|-)\\.(?:n?\\d+|-))\\s*\\)\\s*(.*)");
    static private final Pattern TC_PATTERN = Pattern.compile("(.+?)\\s*\\(\\s*T\\.?C\\.?(?:\\s+|:)(\\d\\.[A-Z]\\.(?:\\d+|-)\\.(?:\\d+|-)\\.(?:\\d+|-)\\s*)\\)\\s*(.*)");
    static private final String PUNCTUATION = ",\\.;:";
    static private final Pattern HYPO_WORD_PATTERN = Pattern.compile("^\\d{7}[a-z]\\d{2}rik\\b|\\b(?:hyphothetical|hyothetical)\\b");

    /**
     * Create a role from its ID and description.
     *
     * @param roleId	a short ID for the role
     * @param roleDesc	the full default description of the role
     */
    public Role(String roleId, String roleDesc) {
        this.roleId = roleId;
        this.description = roleDesc;
        this.checksum = checksum(roleDesc);
    }

    /**
     * @return a checksum used to find this role in a role database
     *
     * @param roleDesc	the role's description; the checksum is the
     * 					same for many small variations in the description
     */
    public static String checksum(String roleDesc) {
        byte[] rText = normalize(roleDesc);
        byte[] digest = DigestUtils.md5(rText);
        String retVal = Base64.encodeBase64String(digest);
        return retVal;
    }

    /**
     * @return the normalized version of a role
     *
     * @param roleDesc
     */
    private static byte[] normalize(String roleDesc) {
        String roleText = roleDesc;
        // Extract the EC and TC numbers.
        String ecNum = null;
        String tcNum = null;
        Matcher m = EC_PATTERN.matcher(roleText);
        if (m.matches()) {
            roleText = join_text(m.group(1), m.group(3));
            ecNum = m.group(2);
        }
        m = TC_PATTERN.matcher(roleText);
        if (m.matches()) {
            roleText = join_text(m.group(1), m.group(3));
            tcNum = m.group(2);
        }
        // Convert to lower case so case doesn't matter.
        roleText = roleText.toLowerCase();
        // Fix spelling mistakes in "hypothetical".
        roleText = HYPO_WORD_PATTERN.matcher(roleText).replaceAll("hypothetical");
        // Remove extra spaces and quotes.
        roleText = roleText.replaceAll("\\r", " ");
        roleText = roleText.trim();
        if (roleText.startsWith("\"")) {
            roleText = roleText.substring(1);
        }
        if (roleText.endsWith("\"")) {
            roleText = roleText.substring(0, roleText.length() - 1);
        }
        roleText = roleText.replaceAll("\\s+", " ");
        // If we have a hypothetical with a number, replace it.
        if (roleText.equals("hypothetical protein") || roleText.isEmpty()) {
            if (ecNum != null) {
                roleText = "putative protein " + ecNum;
            } else if (tcNum != null) {
                roleText = "putative transporter " + tcNum;
            }
        }
        // Now remove the extra spaces and punctuation.
        roleText = roleText.replaceAll("[\\s,.:]{2,}", " ");
        // Convert to UTF-8.
        byte[] retVal = roleText.getBytes(ISO_8859_1);
        return retVal;
    }

    /**
     * Join two text strings.  If there is no punctuation at the start of
     * the second string, a space is inserted.
     *
     * @param string1	first string to use
     * @param string2	string to join to it
     *
     * @return the joined string
     */
    private static String join_text(String string1, String string2) {
        StringBuilder retVal = new StringBuilder(string1.length() +
                string2.length() + 1);
        retVal.append(string1);
        if (string2.length() > 0) {
            char starter = string1.charAt(0);
            if (PUNCTUATION.indexOf(starter) < 0) {
                retVal.append(' ');
            }
            retVal.append(string2);
        }
        return retVal.toString();
    }

    /**
     * This is the sort order for roles.  Two roles with the same description
     * (case-insensitive) will have identical checksums.  If the checksums are
     * equal, the roles are equal, and if the checksums are different, we order
     * by the role description.
     */
    @Override
    public int compareTo(Role o) {
        int retVal = 0;
        if (! this.checksum.equals(o.checksum)) {
            retVal = this.description.compareToIgnoreCase(o.description);
        }
        return retVal;
    }

    /**
     * @return the role's short ID
     */
    public String getRoleId() {
        return roleId;
    }

    /**
     * @return the role's checksum
     */
    public String getChecksum() {
        return checksum;
    }

    /**
     * @return the role's canonical description
     */
    public String getDescription() {
        return description;
    }

    @Override
    public int hashCode() {
        return this.checksum.hashCode();
    }

    /**
     * Note that two roles are equal iff their checksums are equal.
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Role other = (Role) obj;
        if (this.checksum == null) {
            if (other.checksum != null)
                return false;
        } else if (! this.checksum.equals(other.checksum))
            return false;
        return true;
    }

}
