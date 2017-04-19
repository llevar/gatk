package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by shuang on 4/19/17.
 */
public class AlignmentGapBreakerUnitTest {
    @Test
    public void testCompactifyNeighboringSoftClippings() {
        Assert.assertEquals(new Cigar(AlignmentGapBreaker.compactifyNeighboringSoftClippings(TextCigarCodec.decode("1H2S3S4M5D6M7I8M9S10S11H").getCigarElements())),
                TextCigarCodec.decode("1H5S4M5D6M7I8M19S11H"));
    }
}
