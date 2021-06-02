/**
 *
 */
package org.theseed.proteins.kmers.reps;

import java.io.File;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command processes the genomes from a genome source and computes the closest representative genome in one or
 * more RepGen databases.  The output report will list the incoming genome ID and name, and then the ID of the representative
 * in each RepGen database.
 *
 * The positional parameters are the name of the input genome source and the names of the files containing the repgen
 * databases, in order.  Each RepGen database should have a different similarity threshold; otherwise, the column
 * headings will not be unique.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	name of the output file (if not STDOUT)
 * -t	type of input genome source
 *
 * @author Bruce Parrello
 *
 */
public class GtoClassProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoClassProcessor.class);
    /** list of repgen databases */
    private List<RepGenomeDb> repGens;
    /** input genome source */
    private GenomeSource genomes;
    /** list of repgen set names */
    private List<String> repGenNames;

    // COMMAND-LINE OPTIONS

    /** input genome source type */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of input genome source")
    private GenomeSource.Type sourceType;

    /** input genome source file/directory */
    @Argument(index = 0, metaVar = "genomeDir", usage = "input genome source directory", required = true)
    private File genomeDir;

    /** repgen database files */
    @Argument(index = 1, metaVar = "repgen1.ser repgen2.ser ...", usage = "representative-genome databases", required = true)
    private List<File> repGenFiles;

    @Override
    protected void setReporterDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify the repgen files.
        for (File repGenFile : this.repGenFiles) {
            if (! repGenFile.canRead())
                throw new FileNotFoundException("RepGen file " + repGenFile + " is not found or unreadable.");
        }
        // Verify and connect to the genome source.
        if (! this.genomeDir.exists())
            throw new FileNotFoundException("Input genome source " + this.genomeDir + " is not found.");
        log.info("Loading input genomes.");
        this.genomes = this.sourceType.create(this.genomeDir);
        log.info("{} genomes found in {}.", this.genomes.size(), this.genomeDir);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Read in the repGen sets.
        this.buildRepGens(writer);
        // This list will be used to hold the IDs of the representatives found.
        String[] repIds = new String[this.repGens.size()];
        // These sets are used to determine how many groups are found from each repgen set.
        List<Set<String>> groupSets = IntStream.range(0, this.repGenNames.size()).mapToObj(i -> new HashSet<String>(100))
                .collect(Collectors.toList());
        // This will be our progress counter.
        int gCount = 0;
        // These counters track our success rate.
        int foundCount = 0;
        int errorCount = 0;
        // Loop through the genomes.
        for (Genome genome : this.genomes) {
            gCount++;
            log.info("Processing genome {} of {}: {}.", gCount, this.genomes.size(), genome);
            // Find this genome's seed protein.
            String protId = SeqTableProcessor.findSeed(genome);
            if (protId == null) {
                // Here there is no seed protein, so we skip the genome.
                errorCount += repIds.length;
            } else {
                // Get the seed protein kmers.
                Feature seed = genome.getFeature(protId);
                ProteinKmers seedKmers = new ProteinKmers(seed.getProteinTranslation());
                // Find all the representatives.
                for (int i = 0; i < this.repGens.size(); i++) {
                    RepGenomeDb.Representation rep = this.repGens.get(i).findClosest(seedKmers);
                    if (rep.isRepresented()) {
                        // Here we succeeded.
                        repIds[i] = rep.getGenomeId();
                        foundCount++;
                        groupSets.get(i).add(rep.getGenomeId());
                    } else {
                        repIds[i] = genome.getId();
                        errorCount++;
                    }
                }
                // Write this genome's record.
                List<String> columns = new ArrayList<String>(repIds.length + 2);
                columns.add(genome.getId());
                columns.add(genome.getName());
                for (String repId : repIds)
                    columns.add(repId);
                writer.println(StringUtils.join(columns, '\t'));
            }
        }
        log.info("{} genomes processed. {} representatives found; {} failures.", gCount, foundCount, errorCount);
        if (log.isInfoEnabled()) {
            // Describe the number of groups for each repGen set.
            for (int i = 0; i < this.repGenNames.size(); i++)
                log.info("{} groups found using {}.", groupSets.get(i).size(), this.repGenNames.get(i));
        }
    }

    /**
     * Read the repgen sets into memory and create the output header.
     *
     * @param writer	output writer for the report
     * @throws IOException
     */
    private void buildRepGens(PrintWriter writer) throws IOException {
        // We will put the column headings in here.
        StringBuilder headers = new StringBuilder(20 + 12 * this.repGenFiles.size());
        headers.append("genome_id\tgenome_name");
        // This will hold the names of the repgen sets.
        this.repGenNames = new ArrayList<String>(this.repGenFiles.size());
        // This will hold the repgen sets.
        this.repGens = new ArrayList<RepGenomeDb>(this.repGenFiles.size());
        // Build each RepGen set and add its title to the header string.
        for (File repGenFile : this.repGenFiles) {
            log.info("Loading representative genome database from {}.", repGenFile);
            RepGenomeDb repDb = RepGenomeDb.load(repGenFile);
            this.repGens.add(repDb);
            String rgName = String.format("RepGen.%d", repDb.getThreshold());
            this.repGenNames.add(rgName);
            headers.append('\t');
            headers.append(rgName);
        }
        // Write the header line.
        writer.println(headers.toString());
    }

}
