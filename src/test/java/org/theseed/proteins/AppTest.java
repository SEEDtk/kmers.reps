package org.theseed.proteins;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for simple App.
 */
public class AppTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * Test the Role object
     */
    public void testRole() {
        // We need to verify that the normalization works.  We create
        // equivalent roles and insure their checksums are equal.
        String rDesc = "(R)-2-hydroxyacid dehydrogenase, similar to L-sulfolactate dehydrogenase (EC 1.1.1.272)";
        Role rObj1 = new Role("2HydrDehySimiLSulf", rDesc);
        assertEquals("Role ID not stored properly.", "2HydrDehySimiLSulf", rObj1.getRoleId());
        Role rObj2 = new Role("2HydrDehySimiLSulf2", "(R)-2-hydroxyacid dehydrogenase, similar to L-sulfolactate dehydrogenase");
        assertEquals("EC number affects checksum.", rObj1.getChecksum(), rObj2.getChecksum());
        assertEquals("Equal checksums is not equal roles.", rObj1, rObj2);
        assertEquals("Equal checksums is not compare-0 roles.", 0, rObj1.compareTo(rObj2));
        assertEquals("Equal roles have different hashcodes.", rObj1.hashCode(), rObj2.hashCode());
        Role rObj3 = new Role("2HydrDehySimiLSulf3", "(r)-2-hydroxyacid dehydrogenase, similar to L-sulfolactate dehydrogenase");
        assertEquals("Role checksums are case sensitive.", rObj1, rObj3);
        Role rObj4 = new Role("2HydrDehySimiLSulf4", "(R)-2-hydroxyacid dehydrogenase (EC 1.1.1.272), similar to L-sulfolactate dehydrogenase");
        assertEquals("Role checksums affected by EC position.", rObj1, rObj4);
        Role rObj5 = new Role("2HydrDehySimiLSulf5", "(R)-2-hydroxyacid dehydrogenase similar to L-sulfolactate dehydrogenase");
        assertEquals("Role checksum affected by comma.", rObj1, rObj5);
    }
}
