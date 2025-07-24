package org.theseed.proteins.kmers.reps;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

public class CopyTest {

    @Test
    public void testRepDbCopy() throws IOException {
        File rep10file = new File("data", "rep50.db");
        RepGenomeDb repDb1 = RepGenomeDb.load(rep10file);
        RepGenomeDb repDb2 = new RepGenomeDb(100, "Phenylalanyl-tRNA synthetase alpha chain");
        int size = repDb1.size();
        assertThat(size, greaterThan(0));
        BaseGenomeProcessor.copyReps(repDb1, repDb2);
        assertThat(repDb2.size(), equalTo(size));
        for (RepGenome rep1 : repDb1) {
            String rep1Id = rep1.getGenomeId();
            RepGenome rep2 = repDb2.get(rep1Id);
            assertThat(rep1Id, rep2, not(nullValue()));
            assertThat(rep1Id, rep2.getGenomeId(), equalTo(rep1Id));
        }
    }

}
