package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVariantDiscoveryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by shuang on 4/18/17.
 */
public class AlignmentGapBreaker {

    // TODO: 1/18/17 handle mapping quality and mismatches of the new alignment better rather than artificial ones. this would ultimately require a re-alignment done right.
    /**
     * Break a gapped alignment into multiple alignment regions, based on input sensitivity: i.e. if the gap(s) in the input alignment
     * is equal to or larger than the size specified by {@code sensitivity}, two alignment regions will be generated.
     *
     * To fulfill this functionality, needs to accomplish three key tasks correctly:
     * <ul>
     *     <li>generate appropriate cigar,</li>
     *     <li>infer reference coordinates,</li>
     *     <li>infer contig coordinates</li>
     * </ul>
     * As an example, an alignment with CIGAR "397S118M2D26M6I50M7I26M1I8M13D72M398S", when {@code sensitivity} is set to "1", should be broken into 5 alignments
     * <ul>
     *     <li>"397S118M594S"</li>
     *     <li>"515S26M568S"</li>
     *     <li>"547S50M512S"</li>
     *     <li>"631S8M470S"</li>
     *     <li>"639S72M398S"</li>
     * </ul>
     * On the other hand, an alignment with CIGAR "10M10D10M60I10M10I10M50D10M", when {@code sensitivity} is set to "50", should be broken into 3 alignments
     * <ul>
     *     <li>"10M10D10M100S"</li>
     *     <li>"80S10M10I10M10S"</li>
     *     <li>"110S10M"</li>
     * </ul>
     * And when an alignment has hard clipping adjacent to soft clippings, e.g. "1H2S3M5I10M20D6M7S8H", it should be broken into alignments with CIGAR resembling the original CIGAR as much as possible:
     * <ul>
     *     <li>"1H2S3M28S8H"</li>
     *     <li>"1H10S10M13S8H"</li>
     *     <li>"1H20S6M7S8H"</li>
     * </ul>
     *
     * @return an iterable of size >= 1. if size==1, the returned iterable contains only the input (i.e. either no gap or hasn't reached sensitivity)
     */
    static Iterable<SAMRecord> breakGappedAlignment(final SAMRecord originalRecord, final int sensitivity) {

        final int contigUnclippedTotalLength = originalRecord.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping( originalRecord.getCigar() );

        final List<CigarElement> cigarElements = AlignmentGapBreaker.checkCigarAndConvertTerminalInsertionToSoftClip(originalRecord.getReadNegativeStrandFlag() ? CigarUtils.invertCigar(originalRecord.getCigar()) : originalRecord.getCigar());
        if (cigarElements.size() == 1) return new ArrayList<>( Collections.singletonList(originalRecord) );

        final List<SAMRecord> result = new ArrayList<>(3); // blunt guess
        final int originalMapQ = originalRecord.getMappingQuality();

        final List<CigarElement> cigarMemoryList = new ArrayList<>();
        final int clippedNBasesFromStart = SVVariantDiscoveryUtils.getNumClippedBases(true, cigarElements);







        return result;
    }


    /**
     * Checks the input CIGAR for assumption that operator 'D' is not immediately adjacent to clipping operators.
     * Then convert the 'I' CigarElement, if it is at either end (terminal) of the input cigar, to a corresponding 'S' operator.
     * Note that we allow CIGAR of the format '10H10S10I10M', but disallows the format if after the conversion the cigar turns into a giant clip,
     * e.g. '10H10S10I10S10H' is not allowed (if allowed, it becomes a giant clip of '10H30S10H' which is non-sense).
     *
     * @return a pair of number of clipped (hard and soft, including the ones from the converted terminal 'I') bases at the front and back of the
     *         input {@code cigarAlongInput5to3Direction}.
     *
     * @throws IllegalArgumentException when the checks as described above fail.
     */
    @VisibleForTesting
    public static List<CigarElement> checkCigarAndConvertTerminalInsertionToSoftClip(final Cigar cigarAlongInput5to3Direction) {

        if (cigarAlongInput5to3Direction.numCigarElements()<2 ) return cigarAlongInput5to3Direction.getCigarElements();

        final List<CigarElement> cigarElements = new ArrayList<>(cigarAlongInput5to3Direction.getCigarElements());
        SVVariantDiscoveryUtils.validateCigar(cigarElements);

        final List<CigarElement> convertedList = convertInsToSoftClipFromOneEnd(cigarElements, true);
        return convertInsToSoftClipFromOneEnd(convertedList, false);
    }

    /**
     * Actually convert terminal 'I' to 'S' and in case there's an 'S' comes before 'I', compactify the two neighboring 'S' operations into one.
     *
     * @return the converted and compactified list of cigar elements
     */
    @VisibleForTesting
    static List<CigarElement> convertInsToSoftClipFromOneEnd(final List<CigarElement> cigarElements,
                                                             final boolean fromStart) {
        final int numHardClippingBasesFromOneEnd = SVVariantDiscoveryUtils.getNumHardClippingBases(fromStart, cigarElements);
        final int numSoftClippingBasesFromOneEnd = SVVariantDiscoveryUtils.getNumSoftClippingBases(fromStart, cigarElements);

        final int indexOfFirstNonClippingOperation;
        if (numHardClippingBasesFromOneEnd==0 && numSoftClippingBasesFromOneEnd==0) { // no clipping
            indexOfFirstNonClippingOperation = fromStart ? 0 : cigarElements.size()-1;
        } else if (numHardClippingBasesFromOneEnd==0 || numSoftClippingBasesFromOneEnd==0) { // one clipping
            indexOfFirstNonClippingOperation = fromStart ? 1 : cigarElements.size()-2;
        } else {
            indexOfFirstNonClippingOperation = fromStart ? 2 : cigarElements.size()-3;
        }

        final CigarElement element = cigarElements.get(indexOfFirstNonClippingOperation);
        if (element.getOperator() == CigarOperator.I) {

            cigarElements.set(indexOfFirstNonClippingOperation, new CigarElement(element.getLength(), CigarOperator.S));

            return compactifyNeighboringSoftClippings(cigarElements);
        } else {
            return cigarElements;
        }
    }

    /**
     * Compactify two neighboring soft clippings, one of which was converted from an insertion operation.
     * @return the compactified list of operations
     * @throws IllegalArgumentException if there's un-handled edge case where two operations neighboring each other have
     *                                  the same operator (other than 'S') but for some reason was not compactified into one
     */
    @VisibleForTesting
    public static List<CigarElement> compactifyNeighboringSoftClippings(final List<CigarElement> cigarElements) {
        final List<CigarElement> result = new ArrayList<>(cigarElements.size());
        for (final CigarElement element : cigarElements) {
            final int idx = result.size()-1;
            if (result.isEmpty() || result.get(idx).getOperator()!=element.getOperator()) {
                result.add(element);
            } else {
                Utils.validateArg(result.get(idx).getOperator()==CigarOperator.S && element.getOperator()==CigarOperator.S,
                        "Seeing new edge case where two neighboring operations are having the same operator: " + cigarElements.toString());
                result.set(idx, new CigarElement(result.get(idx).getLength()+element.getLength(), CigarOperator.S));
            }
        }
        return result;
    }
}
