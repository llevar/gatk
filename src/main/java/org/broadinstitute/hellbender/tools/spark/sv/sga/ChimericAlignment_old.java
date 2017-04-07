package org.broadinstitute.hellbender.tools.spark.sv.sga;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVariantDiscoveryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * See documentation in {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment}.
 */
@DefaultSerializer(ChimericAlignment_old.Serializer.class)
public class ChimericAlignment_old {

    public final AlignmentRegion regionWithLowerCoordOnContig;
    public final AlignmentRegion regionWithHigherCoordOnContig;

    public final StrandSwitch strandSwitch;
    public final boolean isForwardStrandRepresentation;

    public final byte[] assemblerGivenContigSeq;
    public final List<String> insertionMappings;

    public enum StrandSwitch {
        NO_SWITCH, FORWARD_TO_REVERSE, REVERSE_TO_FORWARD;
    }

    // TODO: 12/14/16 Skipping simple translocation events
    /**
     * Parse all alignment records for a single locally-assembled contig and generate chimeric alignments if available.
     * Applies certain filters to skip the input alignment regions that are:
     *     1) if the alignment region's mapping quality is below a certain threshold, it is skipped
     *     2) if the alignment region is too small, it is skipped
     * If the alignment region passes the above two filters and the next alignment region could be treated as potential inserted sequence,
     * note down the mapping & alignment information of that region and skip it
     *
     * @param input     made of ({alignmentRegions}, sequence) of a signalling locally-assembled contig
     *
     * @return          the chimeric alignments of this sequence (empty if the sequence does not have any alignments)
     */
    @VisibleForTesting
    public static List<ChimericAlignment_old> fromSplitAlignments_old(final Tuple2<Iterable<AlignmentRegion>, byte[]> input) {

        final byte[] sequence = input._2;
        final List<AlignmentRegion> alignmentRegionsSortedByContigCoord = StreamSupport.stream(input._1().spliterator(), false).sorted(Comparator.comparing(a -> AlignmentRegion.startOfAlignmentInContig(a.samRecord))).collect(Collectors.toList());
        if (alignmentRegionsSortedByContigCoord.size() < 2) {
            return new ArrayList<>();
        }

        final List<ChimericAlignment_old> results = new ArrayList<>(alignmentRegionsSortedByContigCoord.size() - 1);
        final List<String> insertionAlignmentRegions = new ArrayList<>();

        final Iterator<AlignmentRegion> iterator = alignmentRegionsSortedByContigCoord.iterator();

        // fast forward to the first alignment region with high MapQ
        AlignmentRegion current = iterator.next();
        while (mapQualTooLow(current) && iterator.hasNext()) {
            current = iterator.next();
        }

        // then iterate over the AR's in pair to identify CA's.
        while ( iterator.hasNext() ) {
            final AlignmentRegion next = iterator.next();
            if (firstAlignmentIsTooShort(current, next, SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH)) {
                continue;
            } else if (nextAlignmentMayBeNovelInsertion(current, next, SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH)) {
                if (iterator.hasNext()) {
                    insertionAlignmentRegions.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }

            final ChimericAlignment_old ca = new ChimericAlignment_old(current, next, sequence, insertionAlignmentRegions);

            final boolean involvesRefIntervalSwitch = involvesRefPositionSwitch(ca.regionWithLowerCoordOnContig, ca.regionWithHigherCoordOnContig);
            final boolean isNotSimpleTranslocation = isNotSimpleTranslocation(ca.regionWithLowerCoordOnContig, ca.regionWithHigherCoordOnContig, ca.strandSwitch, involvesRefIntervalSwitch);
            if (isNotSimpleTranslocation)
                results.add(ca);

            current = next;
        }

        return results;
    }

    // TODO: 11/22/16 it might also be suitable to consider the reference context this alignment region is mapped to and not simply apply a hard filter (need to think about how to test)
    static boolean mapQualTooLow(final AlignmentRegion next) {
        return next.samRecord.getMappingQuality() < SVConstants.DiscoveryStepConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
    }

    @VisibleForTesting
    public static boolean firstAlignmentIsTooShort(final AlignmentRegion first, final AlignmentRegion second,
                                                   final Integer minAlignLength) {
        final boolean old = first.referenceInterval.size() - SVVariantDiscoveryUtils.overlapOnContig(first, second) < minAlignLength;
        final boolean update = first.samRecord.getAlignmentEnd()-first.samRecord.getAlignmentStart()+1 - SVVariantDiscoveryUtils.overlapOnContig(first, second) < minAlignLength;
        return old;    }

    /**
     * To implement the idea that for two consecutive alignment regions of a contig, the one with higher reference coordinate might be a novel insertion.
     */
    @VisibleForTesting
    public static boolean nextAlignmentMayBeNovelInsertion(final AlignmentRegion current, final AlignmentRegion next,
                                                           final Integer minAlignLength) {
        final  SimpleInterval a = new SimpleInterval(current.samRecord);
        final SimpleInterval b = new SimpleInterval(next.samRecord);

        return mapQualTooLow(next) ||                               // inserted sequence might have low mapping quality
                firstAlignmentIsTooShort(next, current, minAlignLength) ||     // inserted sequence might be very small
                a.contains(b) ||   // one might completely contain the other
                b.contains(a);
    }

    /**
     * Construct a new ChimericAlignment_old from two AlignmentRegions.
     * Assumes {@code regionWithLowerCoordOnContig} has a lower {@link AlignmentRegion#startInAssembledContig} than {@code regionWithHigherCoordOnContig},
     * and input alignment regions are NOT completely enclosed in the other.
     */
    @VisibleForTesting
    public ChimericAlignment_old(final AlignmentRegion regionWithLowerCoordOnContig, final AlignmentRegion regionWithHigherCoordOnContig,
                                 final byte[] contigSequence, final List<String> insertionMappings) {

        final String assemblyId = regionWithLowerCoordOnContig.assemblyId;
        final String contigId = regionWithLowerCoordOnContig.contigId;
        Utils.validateArg(assemblyId.equals(regionWithHigherCoordOnContig.assemblyId), "two alignment regions used to construct chimeric alignment are not from the same local assembly.");
        Utils.validateArg(contigId.equals(regionWithHigherCoordOnContig.contigId), "two alignment regions used to construct chimeric alignment are not from the same assembled contig.");

        this.regionWithLowerCoordOnContig = regionWithLowerCoordOnContig;
        this.regionWithHigherCoordOnContig = regionWithHigherCoordOnContig;
        final SimpleInterval interval_1 = new SimpleInterval( regionWithLowerCoordOnContig.samRecord );
        final SimpleInterval interval_2 = new SimpleInterval( regionWithHigherCoordOnContig.samRecord );
        Utils.validateArg(!interval_1.contains(interval_2) && !interval_2.contains(interval_1),
                "one alignment region contains the other, which is wrong " + regionWithLowerCoordOnContig.toPackedString() + regionWithHigherCoordOnContig.toPackedString());

        this.strandSwitch = determineStrandSwitch(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
        final boolean involvesRefIntervalSwitch = involvesRefPositionSwitch(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
        this.isForwardStrandRepresentation = isForwardStrandRepresentation(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig, this.strandSwitch, involvesRefIntervalSwitch);

        this.assemblerGivenContigSeq = contigSequence;

        this.insertionMappings = insertionMappings;
    }

    @VisibleForTesting
    public static StrandSwitch determineStrandSwitch(final AlignmentRegion first, final AlignmentRegion second) {
        if (first.samRecord.getReadNegativeStrandFlag() == second.samRecord.getReadNegativeStrandFlag()) {
            return StrandSwitch.NO_SWITCH;
        } else {
            return first.samRecord.getReadNegativeStrandFlag() ? StrandSwitch.REVERSE_TO_FORWARD : StrandSwitch.FORWARD_TO_REVERSE;
        }
    }

    /**
     * Determine if the region that maps to a lower coordinate on the contig also maps to a lower coordinate on the reference.
     */
    @VisibleForTesting
    public static boolean involvesRefPositionSwitch(final AlignmentRegion regionWithLowerCoordOnContig,
                                                    final AlignmentRegion regionWithHigherCoordOnContig) {

        return regionWithHigherCoordOnContig.samRecord.getAlignmentStart() < regionWithLowerCoordOnContig.samRecord.getAlignmentStart();
    }

    // TODO: 12/15/16 simple translocations are defined here and at this time as inter-chromosomal ones and
    //                intra-chromosomal translocations that involve strand switch
    //                to get the 2nd case right, we need evidence flanking both sides of the inversion, and that could be difficult
    // TODO: 1/17/17 this is used for filtering out possible translocations that we are not currently handling, but it overkills some insertions where the inserted sequence maps to chromosomes other than that of the flanking regions
    /**
     * Determine if the two regions indicate a inter-chromosomal translocation or intra-chromosomal translocation that
     * DOES NOT involve a strand switch.
     */
    @VisibleForTesting
    public static boolean isNotSimpleTranslocation(final AlignmentRegion regionWithLowerCoordOnContig,
                                                   final AlignmentRegion regionWithHigherCoordOnContig,
                                                   final StrandSwitch strandSwitch,
                                                   final boolean involvesReferenceIntervalSwitch) {
        // logic is: must be the same reference chromosome for it not to be an inter-chromosomal translocation
        //           and when regions are mapped to the same reference chromosome, there cannot be reference position swap
        final boolean sameChromosome = regionWithLowerCoordOnContig.samRecord.getContig().equals(regionWithHigherCoordOnContig.samRecord.getContig());
        return sameChromosome
                &&
                (strandSwitch!=StrandSwitch.NO_SWITCH
                        || involvesReferenceIntervalSwitch == regionWithLowerCoordOnContig.samRecord.getReadNegativeStrandFlag());
    }

    /**
     * An SV event could be detected from a contig that seem originate from the forward or reverse strand of the reference,
     * besides the annotation that the alignment flanking regions might flank either side of the two breakpoints.
     */
    @VisibleForTesting
    public static boolean isForwardStrandRepresentation(final AlignmentRegion regionWithLowerCoordOnContig,
                                                        final AlignmentRegion regionWithHigherCoordOnContig,
                                                        final StrandSwitch strandSwitch,
                                                        final boolean involvesReferenceIntervalSwitch) {

        if (strandSwitch == StrandSwitch.NO_SWITCH) {
            return (!regionWithLowerCoordOnContig.samRecord.getReadNegativeStrandFlag()) && (!regionWithHigherCoordOnContig.samRecord.getReadNegativeStrandFlag());
        } else {
            return !involvesReferenceIntervalSwitch;
        }
    }

    public Tuple2<SimpleInterval, SimpleInterval> getCoordSortedReferenceIntervals() {
        if (involvesRefPositionSwitch(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig))
            return new Tuple2<>(new SimpleInterval(regionWithHigherCoordOnContig.samRecord), new SimpleInterval(regionWithLowerCoordOnContig.samRecord));
        else
            return new Tuple2<>(new SimpleInterval(regionWithLowerCoordOnContig.samRecord), new SimpleInterval(regionWithHigherCoordOnContig.samRecord));
    }

    protected ChimericAlignment_old(final Kryo kryo, final Input input) {

        this.regionWithLowerCoordOnContig = kryo.readObject(input, AlignmentRegion.class);
        this.regionWithHigherCoordOnContig = kryo.readObject(input, AlignmentRegion.class);

        this.strandSwitch = StrandSwitch.values()[input.readInt()];
        this.isForwardStrandRepresentation = input.readBoolean();

        this.assemblerGivenContigSeq = kryo.readObject(input, byte[].class);

        @SuppressWarnings("unchecked") // trick to limit scope of unchecked cast warnings, see StackOverflow 509076
        final List<String> tempForWarning = (ArrayList<String>) kryo.readObject(input, ArrayList.class);
        this.insertionMappings = tempForWarning;
    }

    @Override
    public String toString() {
        return regionWithLowerCoordOnContig.toString() +
                "\t" +
                regionWithHigherCoordOnContig.toString() +
                "\t" +
                new String(assemblerGivenContigSeq);
    }

    protected void serialize(final Kryo kryo, final Output output) {

        kryo.writeObject(output, regionWithLowerCoordOnContig);
        kryo.writeObject(output, regionWithHigherCoordOnContig);

        output.writeInt(strandSwitch.ordinal());
        output.writeBoolean(isForwardStrandRepresentation);

        kryo.writeObject(output, assemblerGivenContigSeq);

        kryo.writeObject(output, insertionMappings);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ChimericAlignment_old> {
        @Override
        public void write(final Kryo kryo, final Output output, final ChimericAlignment_old chimericAlignmentOld) {
            chimericAlignmentOld.serialize(kryo, output);
        }

        @Override
        public ChimericAlignment_old read(final Kryo kryo, final Input input, final Class<ChimericAlignment_old> klass ) {
            return new ChimericAlignment_old(kryo, input);
        }
    }
}
