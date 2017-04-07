package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.sga.AlignmentRegion;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.stream.Collectors;

/**
 * This tool takes a SAM file containing the alignments of assembled contigs or long reads to the reference
 * and searches it for split alignments indicating the presence of structural variations. To do so the tool parses
 * primary and supplementary alignments; secondary alignments are ignored. To be considered valid evidence of an SV,
 * two alignments from the same contig must have mapping quality 60, and both alignments must have length greater than
 * or equal to minAlignmentLength.
 */
@CommandLineProgramProperties(summary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        oneLineSummary="Parse a SAM file containing contigs or long reads aligned to the reference, and call SVs",
        programGroup = SparkProgramGroup.class)
public final class DiscoverVariantsFromContigAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(DiscoverVariantsFromContigAlignmentsSpark.class);


    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.MAPPED;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final SAMFileHeader header = getHeaderForReads();
        final JavaRDD<Iterable<SAMRecord>> parsedContigAlignments
                = getReads().groupBy(GATKRead::getName).map(Tuple2::_2).map(r -> breakGappedAlignments(r, header));

        makeSenseAndWrite(parsedContigAlignments, discoverStageArgs.fastaReference, ctx.broadcast(getReference()), getAuthenticatedGCSOptions(),
                vcfOutputFileName, localLogger);
//
//        final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsIterable
//                = getReads().groupBy(GATKRead::getName).map(Tuple2::_2).mapToPair(AssemblyAlignmentParser::breakGappedAlignments);
//
//        makeSenseAndWrite_old(alignmentRegionsIterable, discoverStageArgs.fastaReference, ctx.broadcast(getReference()), getAuthenticatedGCSOptions(),
//                vcfOutputFileName, localLogger);
    }

    public static void makeSenseAndWrite(final JavaRDD<Iterable<SAMRecord>> alignmentRegionsIterable,
                                         final String fastaReference, final Broadcast<ReferenceMultiSource> broadcastReference,
                                         final PipelineOptions pipelineOptions, String vcfFileName,
                                         final Logger toolLogger) {

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments(alignmentRegionsIterable, toolLogger)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(pipelineOptions, vcfFileName, fastaReference, variants, toolLogger);
    }

    /**
     * Converts an iterable of a {@link GATKRead}'s of a contig, particularly the alignment and sequence information in it
     * to the custom {@link AlignmentRegion} format.
     */
    @VisibleForTesting
    static Iterable<SAMRecord> breakGappedAlignments(final Iterable<GATKRead> reads, final SAMFileHeader header) {

        Utils.validateArg(reads.iterator().hasNext(), "input collection of GATK reads is empty");

        return Utils.stream(reads)
                .filter(r -> !r.isSecondaryAlignment())
                .map(r -> AlignmentGapBreaker.breakGappedAlignment(r.convertToSAMRecord(header), SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY))
                .flatMap(Utils::stream).collect(Collectors.toList());
    }

    private static void makeSenseAndWrite_old(final JavaPairRDD<Iterable<AlignmentRegion>, byte[]> alignmentRegionsIterable,
                                              final String fastaReference, final Broadcast<ReferenceMultiSource> broadcastReference,
                                              final PipelineOptions pipelineOptions, String vcfFileName,
                                              final Logger toolLogger) {

        final JavaRDD<VariantContext> variants
                = SVVariantConsensusDiscovery.discoverNovelAdjacencyFromChimericAlignments_old(alignmentRegionsIterable, toolLogger)
                .map(tuple2 -> SVVariantConsensusDiscovery.discoverVariantsFromConsensus(tuple2, broadcastReference));

        SVVCFWriter.writeVCF(pipelineOptions, vcfFileName, fastaReference, variants, toolLogger);
    }
}
