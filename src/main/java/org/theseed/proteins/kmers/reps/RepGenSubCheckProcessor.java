package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.regex.Pattern;

import org.slf4j.LoggerFactory;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This command finds all the repgen databases in a P3Eval directory and verifies each is a subset of the next
 * larger one. It also checks the good-genomes report to insure that all the representative genomes in it make
 * sense.
 * 
 * The positional parameters are the name of the P3Eval directory containing the repgen databases and the
 * name of the good-genomes report file. The repgen databases are expected to be named "rep<threshold>.ser".
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -D   output directory for reports (default: current directory)
 * 
 */
public class RepGenSubCheckProcessor extends BaseMultiReportProcessor {

    // FIELDS
    /** logging facility */
    protected static final Logger log = LoggerFactory.getLogger(RepGenSubCheckProcessor.class);
    /** array of repgen database files */
    private File[] repgenFiles;
    /** array of repgen databases */
    private RepGenomeDb[] repDBs;
    /** filename pattern for repgen databases */
    private static final Pattern REPGEN_PATTERN = Pattern.compile("rep\\d+.ser");
    /** file filter for repgen databases */
    private static final FileFilter REPGEN_FILTER = new FileFilter() {
        @Override
        public boolean accept(File pathname) {
            boolean retVal = pathname.isFile();
            if (retVal) {
                String name = pathname.getName();
                retVal = REPGEN_PATTERN.matcher(name).matches();
            }
            return retVal;
        }
    };

    // COMMAND-LINE OPTIONS

    /** directory containing the repgen files */
    @Argument(index = 0, metaVar = "P3EvalDir", usage = "P3Eval directory containing repgen databases", required = true)
    private File p3EvalDir;

    /** file containing the good-genomes report */
    @Argument(index = 1, metaVar = "bvbrc.good.tbl", usage = "good-genomes report file", required = true)
    private File goodGenomesFile;

    @Override
    protected File setDefaultOutputDir(File curDir) {
        return curDir;
    }

    @Override
    protected void setMultiReportDefaults() {
    }

    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        // Verify the P3Eval directory.
        if (!this.p3EvalDir.isDirectory())
            throw new FileNotFoundException("The specified P3Eval directory " + this.p3EvalDir.getAbsolutePath() + " does not exist or is not a directory.");
        // Verify the good-genomes file.
        if (! this.goodGenomesFile.canRead())
            throw new IOException("The good-genomes report file " + this.goodGenomesFile.getAbsolutePath() + " is not readable.");
        // Get the repgen files.
        this.repgenFiles = this.p3EvalDir.listFiles(REPGEN_FILTER);
        if (this.repgenFiles == null || this.repgenFiles.length == 0)
            throw new IOException("The specified P3Eval directory " + this.p3EvalDir.getAbsolutePath() + " does not contain any repgen databases.");
        // Verify that the files are readable.
        for (File repgenFile : this.repgenFiles) {
            if (!repgenFile.canRead())
                throw new IOException("The repgen database file " + repgenFile.getAbsolutePath() + " is not readable.");
        }
    }

    @Override
    protected void runMultiReports() throws Exception {
        // Now we need to load the repgen databases so we can compare them.
        log.info("Loading repgen databases from {}.", this.p3EvalDir.getAbsolutePath());
        this.repDBs = new RepGenomeDb[this.repgenFiles.length];
        for (int i = 0; i < this.repgenFiles.length; i++) {
            log.info("Loading repgen database {}.", this.repgenFiles[i].getName());
            this.repDBs[i] = RepGenomeDb.load(this.repgenFiles[i]);
        }
        // Now we need to sort the databases by similarity threshold.
        Arrays.sort(this.repDBs, Comparator.comparingInt(x -> x.getThreshold()));
        // First, we do the subset check.
        log.info("Writing subset report.");
        try (PrintWriter writer = this.openReport("repgen.subset.errors.tbl")) {
            // Write the report header.
            writer.println("rep_id\trep_name\trep_db\tmissing_in_db");
            // The procedure is now simple. For each repgen database, we verify all its representative genomes
            // are present in the next database in the array.
            int totalCount = 0;
            int missingCount = 0;        
            for (int i = 0; i < this.repDBs.length - 1; i++) {
                RepGenomeDb db1 = this.repDBs[i];
                RepGenomeDb db2 = this.repDBs[i + 1];
                // Get the names the databases.
                String db1Name = "rep" + db1.getThreshold();
                String db2Name = "rep" + db2.getThreshold();
                log.info("Checking subset relationship between {} and {}.", db1Name, db2Name);
                for (RepGenome rep : db1) {
                    totalCount++;
                    // If the genome is not in the next database, write it to the report.
                    String repId = rep.getGenomeId();
                    if (! db2.contains(repId)) {
                        writer.println(repId + "\t" + rep.getName() + "\t" + db1Name + "\t" + db2Name);
                        missingCount++;
                    }
                }
            }
            log.info("Checked {} genomes in {} repgen databases, found {} missing.", totalCount, this.repDBs.length, missingCount);
        }
        // Now we check the good-genomes report. This requires opening the good-genomes file as well as the output
        // report file.
        log.info("Checking good-genomes report in {}.", this.goodGenomesFile.getAbsolutePath());
        try (PrintWriter writer = this.openReport("repgen.good.errors.tbl");
                TabbedLineReader goodGenomesStream = new TabbedLineReader(this.goodGenomesFile)) {
            // Get the input column indices.
            int idColIdx = goodGenomesStream.findField("genome_id");
            int nameColIdx = goodGenomesStream.findField("name");
            // Now we need the representative-genome column indices.
            int[] repCols = new int[this.repDBs.length];
            for (int i = 0; i < this.repDBs.length; i++)
                repCols[i] = goodGenomesStream.findField("rep" + this.repDBs[i].getThreshold());
            // Initialize the counters.
            int totalCount = 0;
            int invalidRepCount = 0;
            int notSelfCount = 0;
            // Write the report header. We recognize two types of errors: a representative-genome ID that
            // is not in the repgen database, and a representative-genome that does not represent itself.
            writer.println("genome_id\tname\trep_level\tbad_rep_id\terror_type");
            for (TabbedLineReader.Line line : goodGenomesStream) {
                String genomeId = line.get(idColIdx);
                String genomeName = line.get(nameColIdx);
                totalCount++;
                // Check each representative genome.
                for (int i = 0; i < repCols.length; i++) {
                    int repCol = repCols[i];
                    String repId = line.get(repCol);
                    RepGenomeDb repDb = this.repDBs[i];
                    if (! repDb.contains(repId)) {
                        // The representative genome is not in the database.
                        writer.println(genomeId + "\t" + genomeName + "\t" + repDb.getThreshold() + "\t" + repId + "\tNot a Representative");
                        invalidRepCount++;
                    } else if (repDb.contains(genomeId)) {
                        // The representative genome is in the database, but does it represent itself?
                        if (! genomeId.equals(repId)) {
                            // The representative genome is in the database, but it does not represent itself.
                            writer.println(genomeId + "\t" + genomeName + "\t" + repDb.getThreshold() + "\t" + repId + "\tNot Self-Representing");
                            notSelfCount++;
                        }
                    }
                }
            }
            log.info("Checked {} genomes in good-genomes report, found {} invalid representatives and {} not self-representing.", 
                     totalCount, invalidRepCount, notSelfCount);
        }
    }

}
