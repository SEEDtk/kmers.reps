package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This commands reads the good-genomes output from GenomeProcessor and outputs the genomes that do not
 * have hierarchical membership. The repgen sets are designed so that each large-distance representative
 * is also a small-distance representative. Thus, the repgen sets are arranged in a hierarchy like
 * taxonomic groupings. A genome has hierarchical membership if each of its smaller-distance
 * representatives is in the set represented by its larger-distance representatives. The output will be
 * a list of each non-hierarchical genome in "hierarchy.errors.tbl" and an indented hierarchy list in
 * "hierarchy.tree.tbl".
 * 
 * This requires two passes through the good-genomes report: one for creating the hierarchy and one
 * for checking it against the repgen assignments.
 * 
 * The positional parameter is the name of the good-genomes report.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -D   output directory for reports (defaults to current directory)
 */
public class HierarchyCheckProcessor extends BaseMultiReportProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(HierarchyCheckProcessor.class);
    /** map of representative genomes to parent representatives */
    private Map<String, String[]> hierarchyMap;
    /** list of repgen set names, in column order */
    private String[] setNames;
    /** index of the genome ID column in the input file */
    private int idColIdx;
    /** index of the genome name column in the input file */
    private int nameColIdx;
    /** index of the first repgen ID column in the input file */
    private int repCol1Idx;
    /** index of the last repgen ID column in the input file */
    private int lastRepIdx;
    /** descendant maps for each repgen level (no map for the last level) */
    private List<Map<String, Set<String>>> childMapList;
    /** map of repgen genome IDs to names */
    private Map<String, String> nameMap;
    /** map of repgen genome IDs to parents for each level (no map for the first level) */
    private List<Map<String, Set<String>>> parentMapList;
    /** indentation spacer for tree view */
    private static final String SPACER = "    ";

    // COMMAND-LINE OPTIONS

    /** input file name */
    @Argument(index = 0, metaVar = "bvbrc.good.tbl", usage = "good-genomes file containing repgen data", required = true)
    private File inFile;

    @Override
    protected File setDefaultOutputDir(File curDir) {
        return curDir;
    }

    @Override
    protected void setMultiReportDefaults() {
    }

    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        // Create the big data structures.
        this.hierarchyMap = new HashMap<String, String[]>();
        this.childMapList = new ArrayList<Map<String, Set<String>>>();
        this.parentMapList = new ArrayList<Map<String, Set<String>>>();
        this.nameMap = new HashMap<String, String>();
    }

    @Override
    protected void runMultiReports() throws Exception {
        // Make the first pass, where we build the hierarchy map and the tree list.
        log.info("Reading {} to build hierarchy.", this.inFile);
        try (TabbedLineReader inStream = new TabbedLineReader(this.inFile)) {
            // Find the genome-identifying columns.
            this.idColIdx = inStream.findField("genome_id");
            this.nameColIdx = inStream.findField("name");
            // Now find the repgen ID fields.
            String[] colNames = inStream.getLabels();
            int i = 0;
            final int n = colNames.length;
            while (i < n && ! StringUtils.startsWith(colNames[i], "rep")) i++;
            if (i >= n)
                throw new IOException("Missing repgen ID columns in " + this.inFile);
            // Save the repgen labels.
            this.repCol1Idx = i;
            this.setNames = Arrays.copyOfRange(colNames, i, n);
            this.lastRepIdx = n - 1;
            // Initialize the tree lists. Note that the parent map list skips the
            // first level, while the child map list has no last level.
            final int nParentLevels = this.setNames.length - 1;
            this.parentMapList.add(null);
            for (i = 0; i < nParentLevels; i++) {
                this.childMapList.add(new HashMap<String, Set<String>>());
                this.parentMapList.add(new HashMap<String, Set<String>>());
            }            
            // Loop through the input file.
            int inCount = 0;
            for (TabbedLineReader.Line line : inStream) {
                inCount++;
                String genomeId = line.get(this.idColIdx);
                if (genomeId.equals(line.get(this.lastRepIdx))) {
                    // Here we have a representative. Save its ancestor list.
                    String[] ancestors = this.getAncestors(line);
                    this.hierarchyMap.put(genomeId, ancestors);
                    // Save its name.
                    this.nameMap.put(genomeId, line.get(this.nameColIdx));
                    // Now update it in the tree map. For each level in the ancestor array,
                    // we put the genome ID in the parent's child list.
                    String child = genomeId;
                    for (i = nParentLevels - 1; i >= 0; i--) {
                        String parent = ancestors[i];
                        // Get the child set for the parent at this level.
                        Set<String> childSet = this.childMapList.get(i).computeIfAbsent(parent, x -> new TreeSet<String>());
                        // Put the child into it.
                        childSet.add(child);
                        // Get the parent set for the child at the next level.
                        final int iChild = i + 1;
                        Set<String> parentSet = this.parentMapList.get(iChild).computeIfAbsent(child, x -> new TreeSet<String>());
                        // Put the parent into it.
                        parentSet.add(parent);
                        // Set up for the next level. Note that a repgen can have multiple parents this way.
                        child = parent;
                    }
                }
            }
            log.info("{} genomes read, {} representatives found.", inCount, this.hierarchyMap.size());
        }
        // The logic for the hierarchy check is simple: each genome must have the same ancestor list
        // as its representative. We check them one at a time so we can output the smallest-distance
        // one that doesn't match.
        log.info("Re-reading {} to create the bad-hierarchy report.", this.inFile);
        try (TabbedLineReader inStream = new TabbedLineReader(this.inFile); PrintWriter writer = this.openReport("hierarchy.errors.tbl")) {
            writer.println("genome_id\tgenome_name\tset_name\texpected\tactual");
            int inCount = 0;
            int errCount = 0;
            for (TabbedLineReader.Line line : inStream) {
                inCount++;
                String genomeId = line.get(this.idColIdx);
                String lastRep = line.get(this.lastRepIdx);
                String[] ancestors = this.getAncestors(line);
                String[] repAncestors = this.hierarchyMap.get(lastRep);
                if (repAncestors == null)
                    throw new IllegalStateException("Did not find expected ancestor array for " + lastRep + " while checking " + genomeId + ".");
                // Loop backward from the end until we find a mismatch.
                int i = ancestors.length - 1;
                while (i >= 0 && ancestors[i].equals(repAncestors[i])) i--;
                if (i >= 0) {
                    // Here we found a mismatch.
                    errCount++;
                    String errCol = this.setNames[i];
                    writer.println(genomeId + "\t" + line.get(this.nameColIdx) + "\t" + errCol + "\t"
                            + repAncestors[i] + "\t" + ancestors[i]);
                }
            }
            log.info("{} errors found in {} genomes.", errCount, inCount);
        }
        // Next we print the tree. This involves recursion, since we do it depth-first.
        log.info("Writing repgen hierarchy tree.");
        try (PrintWriter writer = this.openReport("hierarchy.tree.txt")) {
            // Now we loop through the high-level genomes.
            Map<String, Set<String>> rootMap = this.childMapList.get(0);
            for (String rootId : rootMap.keySet())
                this.processParent(writer, rootId, 0);
        }
        // Finally, we need a report on the repgens that have ambiguous parents. We process each of the child levels.
        log.info("Writing multiple-parent report.");
        try (PrintWriter writer = this.openReport("hierarchy.multiParent.tbl")) {
            // Since each child will have multiple parents, we will only print the child information on its first line
            // to make the report easier to parse.
            writer.println("rep_id\trep_name\tlevel\tparent_id\tparent_name");
            // Loop through the parent-map list. Note that there is no parent map at
            // the first level.
            int mpCount = 0;
            for (int i = 1; i < this.parentMapList.size(); i++) {
                String levelName = this.setNames[i];
                Map<String, Set<String>> parentMap = this.parentMapList.get(i);
                for (var parentEntry : parentMap.entrySet()) {
                    Set<String> parentSet = parentEntry.getValue();
                    // We only care if there are multiple parents.
                    if (parentSet.size() > 1) {
                        mpCount++;
                        // Print the first parent.
                        String childId = parentEntry.getKey();
                        String childName = this.nameMap.get(childId);
                        Iterator<String> iter = parentSet.iterator();
                        String parentId = iter.next();
                        String parentName = this.nameMap.get(parentId);
                        writer.println(childId + "\t" + childName + "\t" + levelName + "\t" + parentId + "\t" + parentName);
                        // Print the remaining parents.
                        while (iter.hasNext()) {
                            parentId = iter.next();
                            parentName = this.nameMap.get(parentId);
                            writer.println("\t\t\t" + parentId + "\t" + parentName);
                        }
                    }
                }
            }
            log.info("{} repgens found with multiple parents.", mpCount);
        }
    }

    /**
     * Extract the ancestor list from an input line.
     * 
     * @param line      source input line
     * 
     * @return an array of the ancestor IDs, parallel to the setNames array
     */
    private String[] getAncestors(TabbedLineReader.Line line) {
        String[] fields = line.getFields();
        String[] retVal = Arrays.copyOfRange(fields, this.repCol1Idx, fields.length);
        return retVal;
    }

    /**
     * Display the data for a genome and its children.
     * 
     * @param writer        output writer for tree
     * @param parentId      ID of the parent node to output
     * @param parentLevel   level index of the parent
     */
    private void processParent(PrintWriter writer, String parentId, int parentLevel) {
        // Write the parent line.
        String parentName = this.nameMap.get(parentId);
        writer.println(StringUtils.repeat(SPACER, parentLevel) + parentId + " " + parentName);
        // Are there children?
        if (parentLevel < this.childMapList.size()) {
            // Get the child list.
            Map<String, Set<String>> treeMap = this.childMapList.get(parentLevel);
            Set<String> childSet = treeMap.get(parentId);
            if (childSet == null)
                log.warn("No child set found for {}.", parentId);
            else {
                final int childLevel = parentLevel + 1;
                for (String childId : childSet)
                    this.processParent(writer, childId, childLevel);
            }
        }
    }

}
