package org.broadinstitute.hellbender.tools.spark.sv.sga;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentGapBreaker;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Implements a parser for parsing alignments of locally assembled contigs.
 */
public class AssemblyAlignmentParser implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Converts an iterable of a {@link GATKRead}'s of a contig, particularly the alignment and sequence information in it
     * to the custom {@link AlignmentRegion} format.
     */
    @VisibleForTesting
    public static Tuple2<Iterable<AlignmentRegion>, byte[]> convertToAlignmentRegions(final Iterable<GATKRead> reads) {
        Utils.validateArg(reads.iterator().hasNext(), "input collection of GATK reads is empty");

        final GATKRead primaryAlignment = Utils.stream(reads).filter(r -> !(r.isSecondaryAlignment() || r.isSupplementaryAlignment()))
                .findFirst()
                .orElseThrow(() -> new GATKException("no primary alignment for read " + reads.iterator().next().getName()));

        final byte[] bases = primaryAlignment.getBases();
        if (primaryAlignment.isReverseStrand()) {
            SequenceUtil.reverseComplement(bases);
        }

        final Iterable<AlignmentRegion> alignmentRegionIterable
                = Utils.stream(reads)
                .filter(r -> !r.isSecondaryAlignment())
                .map(AlignmentRegion::new)
                .map(ar -> AlignmentGapBreaker.breakGappedAlignment(ar, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY))
                .flatMap(Utils::stream).collect(Collectors.toList());

        return new Tuple2<>(alignmentRegionIterable, bases);
    }

    /**
     * Filter out "failed" assemblies and turn the alignments of contigs into custom {@link AlignmentRegion} format.
     */
    @VisibleForTesting
    public static List<Tuple2<Iterable<AlignmentRegion>, byte[]>> formatToAlignmentRegions(final Iterable<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseIterable,
                                                                                           final List<String> refNames) {

        return Utils.stream(alignedAssemblyOrExcuseIterable)
                .filter(alignedAssemblyOrExcuse -> alignedAssemblyOrExcuse.getErrorMessage() == null)
                .map(alignedAssembly -> forEachAssemblyNotExcuse(alignedAssembly, refNames))
                .flatMap(Utils::stream)                         // size == total # of contigs' from all successful assemblies
                .filter(pair -> pair._1.iterator().hasNext())   // filter out unmapped and contigs without primary alignments
                .collect(Collectors.toList());
    }

    /**
     * Work on "successful" assembly and turn its contigs' alignments to custom {@link AlignmentRegion} format.
     */
    @VisibleForTesting
    public static Iterable<Tuple2<Iterable<AlignmentRegion>, byte[]>> forEachAssemblyNotExcuse(final AlignedAssemblyOrExcuse alignedAssembly,
                                                                                        final List<String> refNames) {

        final FermiLiteAssembly assembly = alignedAssembly.getAssembly();

        final String assemblyIdString = AlignedAssemblyOrExcuse.formatAssemblyID(alignedAssembly.getAssemblyId());

        final List<List<BwaMemAlignment>> allAlignments = alignedAssembly.getContigAlignments();

        final List<Tuple2<Iterable<AlignmentRegion>, byte[]>> result = new ArrayList<>(allAlignments.size());

        IntStream.range(0, assembly.getNContigs())
                .forEach( contigIdx -> {
                    final byte[] contigSequence = assembly.getContig(contigIdx).getSequence();
                    final Iterable<AlignmentRegion> arOfAContig = convertToAlignmentRegions(assemblyIdString, AlignedAssemblyOrExcuse.formatContigID(contigIdx), refNames, allAlignments.get(contigIdx), contigSequence.length);
                    result.add(new Tuple2<>(arOfAContig, contigSequence));
                } );

        return result;
    }

    /**
     * Converts alignment records of the contig pointed to by {@code contigIdx} in a {@link FermiLiteAssembly} to custom {@link AlignmentRegion} format.
     */
    @VisibleForTesting
    static Iterable<AlignmentRegion> convertToAlignmentRegions(final String assemblyIdString, final String contigIdString, final List<String> refNames,
                                                               final List<BwaMemAlignment> contigAlignments, final int unClippedContigLength) {

        return contigAlignments.stream()
                .filter( bwaMemAlignment ->  bwaMemAlignment.getRefId() >= 0 && SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwaMemAlignment.getSamFlag())) // mapped and not XA (i.e. not secondary)
                .map(bwaMemAlignment -> new AlignmentRegion(assemblyIdString, contigIdString, unClippedContigLength, bwaMemAlignment, refNames))
                .map(ar -> AlignmentGapBreaker.breakGappedAlignment(ar, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY))
                .flatMap(Utils::stream).collect(Collectors.toList());
    }
}
