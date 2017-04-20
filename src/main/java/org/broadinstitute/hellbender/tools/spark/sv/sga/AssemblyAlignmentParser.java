package org.broadinstitute.hellbender.tools.spark.sv.sga;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedAssembly;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentGapBreaker;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVariantDiscoveryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Implements a parser for parsing alignments of locally assembled contigs.
 */
class AssemblyAlignmentParser implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Converts an iterable of a {@link GATKRead}'s of a contig, particularly the alignment and sequence information in it
     * to the custom {@link AlignedAssembly.AlignmentInterval} format.
     */
    @VisibleForTesting
    static Tuple2<Iterable<AlignedAssembly.AlignmentInterval>, byte[]> convertToAlignmentRegions(final Iterable<GATKRead> reads, final SAMFileHeader header, final int unClippedContigLength) {

        Utils.validateArg(reads.iterator().hasNext(), "input collection of GATK reads is empty");

        final Iterable<AlignedAssembly.AlignmentInterval> alignmentRegionIterable
                = Utils.stream(reads)
                .filter(r -> !r.isSecondaryAlignment())
                .map(r -> r.convertToSAMRecord(header))
                .map(AlignedAssembly.AlignmentInterval::new)
                .map(ar -> breakGappedAlignment(ar, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, unClippedContigLength))
                .flatMap(Utils::stream).collect(Collectors.toList());

        final GATKRead primaryAlignment = Utils.stream(reads).filter(r -> !(r.isSecondaryAlignment() || r.isSupplementaryAlignment()))
                .findFirst()
                .orElseThrow(() -> new GATKException("no primary alignment for read " + reads.iterator().next().getName()));
        final byte[] bases = primaryAlignment.getBases();
        if (primaryAlignment.isReverseStrand()) {
            SequenceUtil.reverseComplement(bases);
        }

        return new Tuple2<>(alignmentRegionIterable, bases);
    }

    /**
     * See documentation in {@link AlignmentGapBreaker}
     */
    @VisibleForTesting
    static Iterable<AlignedAssembly.AlignmentInterval> breakGappedAlignment(final AlignedAssembly.AlignmentInterval oneRegion, final int sensitivity, int unclippedContigLen) {

        final List<CigarElement> cigarElements = AlignmentGapBreaker.checkCigarAndConvertTerminalInsertionToSoftClip(oneRegion.cigarAlong5to3DirectionOfContig);
        if (cigarElements.size() == 1) return new ArrayList<>( Collections.singletonList(oneRegion) );

        final List<AlignedAssembly.AlignmentInterval> result = new ArrayList<>(3); // blunt guess
        final int originalMapQ = oneRegion.mapQual;

        final List<CigarElement> cigarMemoryList = new ArrayList<>();
        final int clippedNBasesFromStart = SVVariantDiscoveryUtils.getNumClippedBases(true, cigarElements);

        final int hardClippingAtBeginning = cigarElements.get(0).getOperator()==CigarOperator.H ? cigarElements.get(0).getLength() : 0;
        final int hardClippingAtEnd = (cigarElements.get(cigarElements.size()-1).getOperator()== CigarOperator.H)? cigarElements.get(cigarElements.size()-1).getLength() : 0;
        final CigarElement hardClippingAtBeginningMaybeNull = hardClippingAtBeginning==0 ? null : new CigarElement(hardClippingAtBeginning, CigarOperator.H);
        int contigIntervalStart = 1 + clippedNBasesFromStart;
        // we are walking along the contig following the cigar, which indicates that we might be walking backwards on the reference if oneRegion.forwardStrand==false
        int refBoundary1stInTheDirectionOfContig = oneRegion.forwardStrand ? oneRegion.referenceInterval.getStart() : oneRegion.referenceInterval.getEnd();
        for (final CigarElement cigarElement : cigarElements) {
            final CigarOperator op = cigarElement.getOperator();
            final int operatorLen = cigarElement.getLength();
            switch (op) {
                case M: case EQ: case X: case S: case H:
                    cigarMemoryList.add(cigarElement);
                    break;
                case I: case D:
                    if (operatorLen < sensitivity) {
                        cigarMemoryList.add(cigarElement);
                        break;
                    }

                    // collapse cigar memory list into a single cigar for ref & contig interval computation
                    final Cigar memoryCigar = new Cigar(cigarMemoryList);
                    final int effectiveReadLen = memoryCigar.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(memoryCigar) - SVVariantDiscoveryUtils.getNumClippedBases(true, memoryCigar);

                    // task 1: infer reference interval taking into account of strand
                    final SimpleInterval referenceInterval;
                    if (oneRegion.forwardStrand) {
                        referenceInterval = new SimpleInterval(oneRegion.referenceInterval.getContig(),
                                refBoundary1stInTheDirectionOfContig,
                                refBoundary1stInTheDirectionOfContig + (memoryCigar.getReferenceLength()-1));
                    } else {
                        referenceInterval = new SimpleInterval(oneRegion.referenceInterval.getContig(),
                                refBoundary1stInTheDirectionOfContig - (memoryCigar.getReferenceLength()-1), // step backward
                                refBoundary1stInTheDirectionOfContig);
                    }

                    // task 2: infer contig interval
                    final int contigIntervalEnd = contigIntervalStart + effectiveReadLen - 1;

                    // task 3: now add trailing cigar element and create the real cigar for the to-be-returned AR
                    cigarMemoryList.add(new CigarElement(unclippedContigLen-contigIntervalEnd-hardClippingAtEnd, CigarOperator.S));
                    if (hardClippingAtEnd != 0) { // be faithful to hard clipping (as the accompanying bases have been hard-clipped away)
                        cigarMemoryList.add(new CigarElement(hardClippingAtEnd, CigarOperator.H));
                    }
                    final Cigar cigarForNewAlignmentRegion = new Cigar(cigarMemoryList);

                    final AlignedAssembly.AlignmentInterval split = new AlignedAssembly.AlignmentInterval(referenceInterval, contigIntervalStart, contigIntervalEnd, cigarForNewAlignmentRegion, oneRegion.forwardStrand, originalMapQ, SVConstants.DiscoveryStepConstants.ARTIFICIAL_MISMATCH);

                    result.add(split);

                    // update cigar memory
                    cigarMemoryList.clear();
                    if (hardClippingAtBeginningMaybeNull != null) {
                        cigarMemoryList.add(hardClippingAtBeginningMaybeNull); // be faithful about hard clippings
                    }
                    cigarMemoryList.add(new CigarElement(contigIntervalEnd - hardClippingAtBeginning + (op.consumesReadBases() ? operatorLen : 0), CigarOperator.S));

                    // update pointers into reference and contig
                    final int refBoundaryAdvance = op.consumesReadBases() ? memoryCigar.getReferenceLength() : memoryCigar.getReferenceLength() + operatorLen;
                    refBoundary1stInTheDirectionOfContig += oneRegion.forwardStrand ? refBoundaryAdvance : -refBoundaryAdvance;
                    contigIntervalStart += op.consumesReadBases() ? effectiveReadLen + operatorLen : effectiveReadLen;

                    break;
                default:
                    throw new GATKException("Alignment CIGAR contains an unexpected N or P element: " + oneRegion.toString()); // TODO: 1/20/17 still not quite sure if this is quite right, it doesn't blow up on NA12878 WGS, but who knows what happens in the future
            }
        }

        if (result.isEmpty()) {
            return new ArrayList<>(Collections.singletonList(oneRegion));
        }

        final SimpleInterval lastReferenceInterval;
        if (oneRegion.forwardStrand) {
            lastReferenceInterval =  new SimpleInterval(oneRegion.referenceInterval.getContig(), refBoundary1stInTheDirectionOfContig, oneRegion.referenceInterval.getEnd());
        } else {
            lastReferenceInterval =  new SimpleInterval(oneRegion.referenceInterval.getContig(), oneRegion.referenceInterval.getStart(), refBoundary1stInTheDirectionOfContig);
        }

        final Cigar lastForwardStrandCigar = new Cigar(cigarMemoryList);
        int clippedNBasesFromEnd = SVVariantDiscoveryUtils.getNumClippedBases(false, cigarElements);
        result.add(new AlignedAssembly.AlignmentInterval(lastReferenceInterval,
                contigIntervalStart, unclippedContigLen-clippedNBasesFromEnd, lastForwardStrandCigar,
                oneRegion.forwardStrand, originalMapQ, SVConstants.DiscoveryStepConstants.ARTIFICIAL_MISMATCH));


        return result;
    }

//    /**
//     * Filter out "failed" assemblies and turn the alignments of contigs into custom {@link AlignmentRegion} format.
//     */
//    @VisibleForTesting
//    static List<Tuple2<Iterable<AlignedAssembly.AlignmentInterval>, byte[]>> formatToAlignmentRegions(final Iterable<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseIterable,
//                                                                                    final List<String> refNames) {
//
//        return Utils.stream(alignedAssemblyOrExcuseIterable)
//                .filter(alignedAssemblyOrExcuse -> alignedAssemblyOrExcuse.getErrorMessage() == null)
//                .map(alignedAssembly -> forEachAssemblyNotExcuse(alignedAssembly, refNames))
//                .flatMap(Utils::stream)                         // size == total # of contigs' from all successful assemblies
//                .filter(pair -> pair._1.iterator().hasNext())   // filter out unmapped and contigs without primary alignments
//                .collect(Collectors.toList());
//    }
//
//    /**
//     * Work on "successful" assembly and turn its contigs' alignments to custom {@link AlignmentRegion} format.
//     */
//    @VisibleForTesting
//    static Iterable<Tuple2<Iterable<AlignedAssembly.AlignmentInterval>, byte[]>> forEachAssemblyNotExcuse(final AlignedAssemblyOrExcuse alignedAssembly,
//                                                                                        final List<String> refNames) {
//
//        final FermiLiteAssembly assembly = alignedAssembly.getAssembly();
//
//        final String assemblyIdString = AlignedAssemblyOrExcuse.formatAssemblyID(alignedAssembly.getAssemblyId());
//
//        final List<List<BwaMemAlignment>> allAlignments = alignedAssembly.getContigAlignments();
//
//        final List<Tuple2<Iterable<AlignedAssembly.AlignmentInterval>, byte[]>> result = new ArrayList<>(allAlignments.size());
//
//        IntStream.range(0, assembly.getNContigs())
//                .forEach( contigIdx -> {
//                    final byte[] contigSequence = assembly.getContig(contigIdx).getSequence();
//                    final Iterable<AlignedAssembly.AlignmentInterval> arOfAContig = convertToAlignmentRegions(assemblyIdString, AlignedAssemblyOrExcuse.formatContigID(contigIdx), refNames, allAlignments.get(contigIdx), contigSequence.length);
//                    result.add(new Tuple2<>(arOfAContig, contigSequence));
//                } );
//
//        return result;
//    }
//
//    /**
//     * Converts alignment records of the contig pointed to by {@code contigIdx} in a {@link FermiLiteAssembly} to custom {@link AlignmentRegion} format.
//     */
//    @VisibleForTesting
//    private static Iterable<AlignedAssembly.AlignmentInterval> convertToAlignmentRegions(final String assemblyIdString, final String contigIdString, final List<String> refNames,
//                                                                       final List<BwaMemAlignment> contigAlignments, final int unClippedContigLength) {
//
//        return contigAlignments.stream()
//                .filter( bwaMemAlignment ->  bwaMemAlignment.getRefId() >= 0 && SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwaMemAlignment.getSamFlag())) // mapped and not XA (i.e. not secondary)
//                .map(bwaMemAlignment -> BwaMemAlignmentUtils.applyAlignment(assemblyIdString+":"+contigIdString, ))
//                .map(AlignedAssembly.AlignmentInterval::new)
//                /*.map(bwaMemAlignment -> new AlignmentRegion(assemblyIdString, contigIdString, unClippedContigLength, bwaMemAlignment, refNames))*/
//                .map(ar -> breakGappedAlignment(ar, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, unClippedContigLength))
//                .flatMap(Utils::stream).collect(Collectors.toList());
//    }

}
