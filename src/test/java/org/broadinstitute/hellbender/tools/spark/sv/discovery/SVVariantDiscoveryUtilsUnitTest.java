package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.sga.AlignmentRegion;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVariantDiscoveryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;


public class SVVariantDiscoveryUtilsUnitTest {

    @Test
    public void testAlignmentRegionOverlap() throws Exception {

        final AlignmentRegion ar1 = new AlignmentRegion("1","1", new SimpleInterval("1",1,5), TextCigarCodec.decode("5M5H"),true, 60, 0, 1,5);
        final AlignmentRegion ar2 = new AlignmentRegion("1","1", new SimpleInterval("1",10,16), TextCigarCodec.decode("4S6M"),true, 60, 0, 5,10);
        Assert.assertEquals(SVVariantDiscoveryUtils.overlapOnContig(ar1, ar2), 1);

        final AlignmentRegion ar3 = new AlignmentRegion("1","1", new SimpleInterval("1",1,5), TextCigarCodec.decode("5M5H"),true, 60, 0, 1,5);
        final AlignmentRegion ar4 = new AlignmentRegion("1","1", new SimpleInterval("1",11,16), TextCigarCodec.decode("5S5M"),true, 60, 0, 6,10);
        Assert.assertEquals(SVVariantDiscoveryUtils.overlapOnContig(ar3, ar4), 0);
    }

    @Test
    public void testClippingArithmetic() {
        Cigar cigar = TextCigarCodec.decode("100M51S");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51S100M");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("100M51H");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51H100M");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("12H12S101M13S13H");
        Assert.assertEquals(SVVariantDiscoveryUtils.getTotalHardClipping(cigar), 25);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(true, cigar), 24);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumClippedBases(false, cigar), 26);
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_emptyCigarElementList(){
        @SuppressWarnings("unchecked")
        final List<CigarElement> emptyList = Collections.EMPTY_LIST;
        SVVariantDiscoveryUtils.validateCigar(emptyList);
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_deletionNeighboringClipping(){
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10S10D10M").getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10H10S10D10M").getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10M10D10S").getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10M10D10S10H").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_only1NonAlignment(){
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10S").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testCigarChecker_noAlignment(){
        SVVariantDiscoveryUtils.validateCigar(TextCigarCodec.decode("10H10S10I10S10H").getCigarElements());
    }

    @Test
    public void testGetNumClippingBases_hardAndSoftSeparately() {
        List<CigarElement> cigarElements = TextCigarCodec.decode("10H20S30M40D50M60S70H").getCigarElements();
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumHardClippingBases(true, cigarElements), 10);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumHardClippingBases(false, cigarElements), 70);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumSoftClippingBases(true, cigarElements), 20);
        Assert.assertEquals(SVVariantDiscoveryUtils.getNumSoftClippingBases(false, cigarElements), 60);
    }

    @Test
    public void testGetIndexOfFirstNonClippingBase(){
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), true), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), false), 0);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10S10D10M").getCigarElements(), true), 1);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10H10S10D10M").getCigarElements(), true), 2);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S").getCigarElements(), false), 1);
        Assert.assertEquals(SVVariantDiscoveryUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S10H").getCigarElements(), false), 1);
    }
}
