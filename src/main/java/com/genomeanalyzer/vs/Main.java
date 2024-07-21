package com.genomeanalyzer.vs;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

import java.net.URL;

public class Main {
    public static void main(String[] args) throws Exception {
        ProteinSequence sequence1 = getSequenceForId("P06213");
        System.out.println(sequence1);
        ProteinSequence sequence2 = getSequenceForId("P14735");
        System.out.println(sequence2);

        alignPairGlobal(sequence1, sequence2);
    }

    private static void alignPairGlobal(ProteinSequence s1, ProteinSequence s2) {
        SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum100();
        SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(s1, s2,
                Alignments.PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(), matrix);
        System.out.printf("%n%s vs %s%n%s", pair.getQuery().getAccession(), pair.getTarget().getAccession(), pair);
    }

    public static ProteinSequence getSequenceForId(String uniprotId) throws Exception {
        URL uniprotFasta = new URL(String.format("https://rest.uniprot.org/uniprotkb/%s.fasta", uniprotId));
        return FastaReaderHelper.readFastaProteinSequence(uniprotFasta.openStream()).get(uniprotId);
    }
}