package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple5;

public class ChimericAlignmentUnitTest extends BaseTest {

    @Test
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = SVDiscoveryTestDataProvider.LONG_CONTIG1.getBytes();
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36);

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test
    public void testFilterByNextAlignmentMayBeNovelInsertion() throws Exception {
        final AlignedAssembly.AlignmentInterval overlappingRegion1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("19", 48699881, 48700035), 1, 154, TextCigarCodec.decode("47S154M"), false, 60, 0);
        final AlignedAssembly.AlignmentInterval overlappingRegion2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("19", 48700584, 48700669), 117, 201, TextCigarCodec.decode("116H85M"), true, 60, 0);

        Assert.assertTrue(ChimericAlignment.nextAlignmentMayBeNovelInsertion(overlappingRegion1, overlappingRegion2, 50));
    }

    @Test
    public void testBooleanStates_inversion() {

        Tuple5<AlignedAssembly.AlignmentInterval, AlignedAssembly.AlignmentInterval, NovelAdjacencyReferenceLocations, String, String> testData = SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.REVERSE_TO_FORWARD, false, true, true);

        testData = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_leftPlus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.FORWARD_TO_REVERSE, false, true, true);

        testData = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_leftMinus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.FORWARD_TO_REVERSE, true, true, false);

        testData = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_rightPlus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.REVERSE_TO_FORWARD, false, true, true);

        testData = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_rightMinus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.REVERSE_TO_FORWARD, true, true, false);
    }

    @Test
    public void testBooleanStates_simpleInsertionAndDeletion() {

        // simple deletion
        Tuple5<AlignedAssembly.AlignmentInterval, AlignedAssembly.AlignmentInterval, NovelAdjacencyReferenceLocations, String, String> testData = SVDiscoveryTestDataProvider.forSimpleDeletion_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // simple insertion
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // long range substitution
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // simple deletion with homology
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);
    }

    @Test
    public void testBooleanStates_tandemDuplication_simple() {

        // tandem duplication simple contraction
        Tuple5<AlignedAssembly.AlignmentInterval, AlignedAssembly.AlignmentInterval, NovelAdjacencyReferenceLocations, String, String> testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // tandem duplication simple expansion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // tandem duplication simple expansion with novel insertion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);
    }

    @Test
    public void testBooleanStates_tandemDuplication_complex() {

        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        Tuple5<AlignedAssembly.AlignmentInterval, AlignedAssembly.AlignmentInterval, NovelAdjacencyReferenceLocations, String, String> testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // third test: contraction from 3 units to 2 units without pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);
    }

    private static void testBooleanSeries(final AlignedAssembly.AlignmentInterval region1, final AlignedAssembly.AlignmentInterval region2,
                                          final ChimericAlignment.StrandSwitch expectedStrandSwitch,
                                          final boolean expectedRefPositionSwitch,
                                          final boolean expectedIsNotSimpleTranslocation,
                                          final boolean expectedIsForwardStrandRepresentation) {

        Assert.assertEquals(ChimericAlignment.determineStrandSwitch(region1, region2), expectedStrandSwitch);
        Assert.assertEquals(ChimericAlignment.involvesRefPositionSwitch(region1, region2), expectedRefPositionSwitch);
        Assert.assertEquals(ChimericAlignment.isForwardStrandRepresentation(region1, region2, expectedStrandSwitch, expectedRefPositionSwitch), expectedIsForwardStrandRepresentation);
        Assert.assertEquals(ChimericAlignment.isNotSimpleTranslocation(region1, region2, expectedStrandSwitch, expectedRefPositionSwitch), expectedIsNotSimpleTranslocation);
    }
}