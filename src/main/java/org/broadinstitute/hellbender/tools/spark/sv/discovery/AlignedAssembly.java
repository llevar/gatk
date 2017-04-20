package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVariantDiscoveryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by shuang on 4/19/17.
 */
@DefaultSerializer(AlignedAssembly.Serializer.class)
public final class AlignedAssembly {

    final int assemblyId;

    final List<AlignedContig> listOfContigsWithItsAlignmentIntervals;

    /**
     * Each assembled contig should have at least one such accompanying structure, or 0 when it is unmapped.
     */
    @DefaultSerializer(AlignmentInterval.Serializer.class)
    static public final class AlignmentInterval {

        public final SimpleInterval referenceInterval;
        public final int startInAssembledContig;   // 1-based, inclusive
        public final int endInAssembledContig;     // 1-based, inclusive

        public final Cigar cigarAlong5to3DirectionOfContig;

        public final boolean forwardStrand;
        public final int mapQual;
        public final int mismatches;

        public AlignmentInterval(final SAMRecord samRecord) {

            final boolean isMappedReverse = samRecord.getReadNegativeStrandFlag();
            this.referenceInterval = new SimpleInterval(samRecord);
            this.startInAssembledContig = getAlignmentStartInOriginalContig(samRecord);
            this.endInAssembledContig = getAlignmentEndInOriginalContig(samRecord);

            this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(samRecord.getCigar()) : samRecord.getCigar();
            this.forwardStrand = !isMappedReverse;
            this.mapQual = samRecord.getMappingQuality();
            final Integer nMismatches = samRecord.getIntegerAttribute("NM");
            this.mismatches = nMismatches==null ? SVConstants.DiscoveryStepConstants.MISSING_NM : nMismatches;
        }

        public AlignmentInterval(final SimpleInterval referenceInterval, final int startInAssembledContig, final int endInAssembledContig,
                                 final Cigar cigarAlong5to3DirectionOfContig, final boolean forwardStrand, final int mapQual, final int mismatches) {
            this.referenceInterval = referenceInterval;
            this.startInAssembledContig = startInAssembledContig;   // 1-based, inclusive
            this.endInAssembledContig = endInAssembledContig;     // 1-based, inclusive

            this.cigarAlong5to3DirectionOfContig = cigarAlong5to3DirectionOfContig;

            this.forwardStrand = forwardStrand;
            this.mapQual = mapQual;
            this.mismatches = mismatches;
        }


        static int getAlignmentStartInOriginalContig(final SAMRecord samRecord) {
            return SVVariantDiscoveryUtils.getNumClippedBases(!samRecord.getReadNegativeStrandFlag(), samRecord.getCigar()) + 1;
        }

        static int getAlignmentEndInOriginalContig(final SAMRecord samRecord) {
            final Cigar cigar = samRecord.getCigar();
            return cigar.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(cigar) - SVVariantDiscoveryUtils.getNumClippedBases(samRecord.getReadNegativeStrandFlag(), cigar);
        }

        AlignmentInterval(final Kryo kryo, final Input input) {
            final String chr = input.readString();
            final int refStart = input.readInt(),
                    refEnd = input.readInt();
            referenceInterval = new SimpleInterval(chr, refStart, refEnd);
            startInAssembledContig = input.readInt();
            endInAssembledContig = input.readInt();
            cigarAlong5to3DirectionOfContig = TextCigarCodec.decode(input.readString());
            forwardStrand = input.readBoolean();
            mapQual = input.readInt();
            mismatches = input.readInt();
        }

        void serialize(final Kryo kryo, final Output output) {
            output.writeString(referenceInterval.getContig());
            output.writeInt(startInAssembledContig);
            output.writeInt(endInAssembledContig);
            output.writeString(TextCigarCodec.encode(cigarAlong5to3DirectionOfContig));
            output.writeBoolean(forwardStrand);
            output.writeInt(mapQual);
            output.writeInt(mismatches);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignmentInterval> {
            @Override
            public void write( final Kryo kryo, final Output output, final AlignmentInterval alignmentInterval){
                alignmentInterval.serialize(kryo, output);
            }

            @Override
            public AlignmentInterval read(final Kryo kryo, final Input input, final Class<AlignmentInterval> clazz ) {
                return new AlignmentInterval(kryo, input);
            }
        }

        // TODO: 11/27/16 test
        /**
         * @return  A packed String representation of this AlignmentRegion, noticeably the fields are not separated by tab.
         *          Note that the format is NOT the same as that used in {@link #toString()}.
         */
        public String toPackedString(final String assemblyId, final String contigId) {
            return String.join(PACKED_STRING_REP_SEPARATOR,
                    assemblyId, contigId, String.valueOf(startInAssembledContig), String.valueOf(endInAssembledContig),
                    referenceInterval.getContig(), String.valueOf(referenceInterval.getStart()), (forwardStrand ? "+" : "-"),
                    TextCigarCodec.encode(cigarAlong5to3DirectionOfContig), String.valueOf(mapQual), String.valueOf(mismatches));
        }
        public static final String PACKED_STRING_REP_SEPARATOR= "_";

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            AlignmentInterval that = (AlignmentInterval) o;

            if (startInAssembledContig != that.startInAssembledContig) return false;
            if (endInAssembledContig != that.endInAssembledContig) return false;
            if (forwardStrand != that.forwardStrand) return false;
            if (mapQual != that.mapQual) return false;
            if (mismatches != that.mismatches) return false;
            if (!referenceInterval.equals(that.referenceInterval)) return false;
            return cigarAlong5to3DirectionOfContig.equals(that.cigarAlong5to3DirectionOfContig);
        }

        @Override
        public int hashCode() {
            int result = referenceInterval.hashCode();
            result = 31 * result + startInAssembledContig;
            result = 31 * result + endInAssembledContig;
            result = 31 * result + cigarAlong5to3DirectionOfContig.hashCode();
            result = 31 * result + (forwardStrand ? 1 : 0);
            result = 31 * result + mapQual;
            result = 31 * result + mismatches;
            return result;
        }
    }

    @DefaultSerializer(AlignedContig.Serializer.class)
    static final class AlignedContig {

        final byte[] contigSequence;
        final List<AlignmentInterval> alignmentIntervals;

        AlignedContig(final byte[] contigSequence, final Iterable<SAMRecord> samRecordList) {
            this.contigSequence = contigSequence;
            this.alignmentIntervals = Utils.stream(samRecordList).map(AlignmentInterval::new).collect(Collectors.toList());
        }

        AlignedContig(final Kryo kryo, final Input input) {

            final int nBases = input.readInt();
            contigSequence = new byte[nBases];
            for(int b=0; b<nBases; ++b) {
                contigSequence[b] = input.readByte();
            }

            final int nAlignments = input.readInt();
            alignmentIntervals = new ArrayList<>(nAlignments);
            for(int i=0; i<nAlignments; ++i) {
                alignmentIntervals.add(new AlignmentInterval(kryo, input));
            }
        }

        void serialize(final Kryo kryo, final Output output) {

            output.writeInt(contigSequence.length);
            for(final byte base : contigSequence) {
                output.writeByte(base);
            }

            output.writeInt(alignmentIntervals.size());
            for(final AlignmentInterval alignmentInterval : alignmentIntervals) {
                alignmentInterval.serialize(kryo, output);
            }
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedContig> {
            @Override
            public void write( final Kryo kryo, final Output output, final AlignedContig alignedContig){
                alignedContig.serialize(kryo, output);
            }

            @Override
            public AlignedContig read(final Kryo kryo, final Input input, final Class<AlignedContig> clazz ) {
                return new AlignedContig(kryo, input);
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            AlignedContig that = (AlignedContig) o;

            if (!Arrays.equals(contigSequence, that.contigSequence)) return false;
            return alignmentIntervals.equals(that.alignmentIntervals);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(contigSequence);
            result = 31 * result + alignmentIntervals.hashCode();
            return result;
        }

        /**
         * todo Remember this does not break the gaps
         * @param samRecords
         * @param unClippedContigLength
         * @return
         */
        static Iterable<AlignedAssembly.AlignmentInterval> extractAlignmentIntervalsFromSAM(final Iterable<SAMRecord> samRecords,
                                                                                            final int unClippedContigLength) {
            return Utils.stream(samRecords).map(AlignedAssembly.AlignmentInterval::new).collect(Collectors.toList());
        }
    }


    AlignedAssembly(final int assemblyId, final FermiLiteAssembly assembly, final List<Iterable<SAMRecord>> allAlignments) {
        Utils.validate(assembly.getNContigs() == allAlignments.size(),
                "input assembly " + String.valueOf(assemblyId) + " doesn't have matching number of contigs " + assembly.getNContigs() + " and number of alignments clusters (assumption that each contig has a alignment cluster--empty if unmapped--is broken).");

        this.assemblyId = assemblyId;
        final int nContigs = assembly.getNContigs();

        this.listOfContigsWithItsAlignmentIntervals
                = IntStream.range(0, nContigs)
                .mapToObj(contigIdx -> new AlignedContig(assembly.getContig(contigIdx).getSequence(), allAlignments.get(contigIdx)))
                .collect(Collectors.toList());
    }

    private AlignedAssembly(final Kryo kryo, final Input input) {
        this.assemblyId = input.readInt();

        final int nContigs = input.readInt();
        listOfContigsWithItsAlignmentIntervals = new ArrayList<>(nContigs);
        for(int contigIdx = 0; contigIdx < nContigs; ++contigIdx) {
            listOfContigsWithItsAlignmentIntervals.add(new AlignedContig(kryo, input));
        }
    }

    private void serialize(final Kryo kryo, final Output output) {
        output.writeInt(assemblyId);

        output.writeInt(listOfContigsWithItsAlignmentIntervals.size());
        for(final AlignedContig alignedContig : listOfContigsWithItsAlignmentIntervals) {
            alignedContig.serialize(kryo, output);
        }
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedAssembly> {
        @Override
        public void write( final Kryo kryo, final Output output, final AlignedAssembly alignedAssembly){
            alignedAssembly.serialize(kryo, output);
        }

        @Override
        public AlignedAssembly read(final Kryo kryo, final Input input, final Class<AlignedAssembly> clazz ) {
            return new AlignedAssembly(kryo, input);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AlignedAssembly that = (AlignedAssembly) o;

        if (assemblyId != that.assemblyId) return false;
        return listOfContigsWithItsAlignmentIntervals.equals(that.listOfContigsWithItsAlignmentIntervals);
    }

    @Override
    public int hashCode() {
        int result = assemblyId;
        result = 31 * result + listOfContigsWithItsAlignmentIntervals.hashCode();
        return result;
    }
}
