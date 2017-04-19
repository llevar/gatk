package org.broadinstitute.hellbender.tools.spark.sv.sga;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import scala.Tuple2;

import java.util.*;

class ContigAligner {

    private final String indexImageFile;

    ContigAligner(final String indexImageFile) {
        this.indexImageFile = indexImageFile;
    }

    /**
     * Takes a collection of assembled contigs and aligns them to the reference with bwa-mem. Non-canonical
     * (secondary) alignments are filtered out, preserving the primary and supplementary alignments.
     * Within the output list, alignments are sorted first by contig (based on the order in which
     * the contigs were passed in, and then by their start position on the contig).
     *
     * @param assemblyId An identifier for the assembly or set of contigs
     * @param contigsCollection The set of all canonical (primary or supplementary) alignments for the contigs.
     */
    List<AlignmentRegion> alignContigs(final String assemblyId, final ContigsCollection contigsCollection) {
        final List<AlignmentRegion> alignedContigs = new ArrayList<>(contigsCollection.getContents().size());
        final BwaMemIndex index = BwaMemIndexSingleton.getInstance(indexImageFile);
        try ( final BwaMemAligner aligner = new BwaMemAligner(index) ) {
            final List<String> refNames = index.getReferenceContigNames();
            final List<Tuple2<ContigsCollection.ContigID, ContigsCollection.ContigSequence>> contents = contigsCollection.getContents();
            final List<byte[]> seqs = new ArrayList<>(contents.size());
            for ( final Tuple2<ContigsCollection.ContigID, ContigsCollection.ContigSequence> contigInfo : contents ) {
                seqs.add(contigInfo._2().toString().getBytes());
            }
            final Iterator<List<BwaMemAlignment>> alignmentsItr = aligner.alignSeqs(seqs).iterator();
            for ( final Tuple2<ContigsCollection.ContigID, ContigsCollection.ContigSequence> contigInfo : contents ) {
                final String contigId = contigInfo._1.toString();
                final int contigLen = contigInfo._2().toString().length();
                final List<BwaMemAlignment> alignments = alignmentsItr.next();

                // filter out secondary alignments, convert to AlignmentRegion objects and sort by alignment start pos
                alignments.stream()
                        .filter(a -> (a.getSamFlag()&SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue())==0)
                        .filter(a -> (a.getSamFlag()&SAMFlag.READ_UNMAPPED.intValue())==0)
                        .map(a -> new AlignmentRegion(assemblyId, contigId, contigLen, a, refNames))
                        .sorted(Comparator.comparing(a -> a.startInAssembledContig))
                        .forEach(alignedContigs::add);
            }
        }

        return alignedContigs;
    }
}
