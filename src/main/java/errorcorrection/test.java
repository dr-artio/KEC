/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package errorcorrection;

import java.io.File;
import java.io.IOException;
import java.util.StringTokenizer;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.*;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;


/**
 *
 * @author kki8
 */
public class test {
    public static void main(String[] args) throws IOException, InterruptedException {

               String targetSeq = "CACGTTTCTTGTGGCAGCTTAAGTTTGAATGTCATTTCTTCAATGGGACGGA"+
		          "GCGGGTGCGGTTGCTGGAAAGATGCATCTATAACCAAGAGGAGTCCGTGCGCTTCGACAGC"+
			  "GACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACA"+
			  "GCCAGAAGGACCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACCTACTGCAGACACAACTA"+ 
			  "CGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAG";
		DNASequence target = new DNASequence(targetSeq,
				AmbiguityDNACompoundSet.getDNACompoundSet());
 
		String querySeq = "ACGAGTGCGTGTTTTCCCGCCTGGTCCCCAGGCCCCCTTTCCGTCCTCAGGAA"+
			  "GACAGAGGAGGAGCCCCTCGGGCTGCAGGTGGTGGGCGTTGCGGCGGCGGCCGGTTAAGGT"+
			  "TCCCAGTGCCCGCACCCGGCCCACGGGAGCCCCGGACTGGCGGCGTCACTGTCAGTGTCTT"+
			  "CTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTAC"+
			  "TCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTCCTGGACAGATACT"+
			  "TCCATAACCAGGAGGAGAACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGAC"+
			  "GGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACATCCTGGAAGACGAG"+
			  "CGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGAGAGCTTCACCGTGCA"+ 
			  "GCGGCGAGACGCACTCGT";
		DNASequence query = new DNASequence(querySeq,
				AmbiguityDNACompoundSet.getDNACompoundSet());
 
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
 
		SimpleGapPenalty gapP = new SimpleGapPenalty();
		gapP.setOpenPenalty((short)15);
		gapP.setExtensionPenalty((short)6);
 
		SequencePair<DNASequence, NucleotideCompound> psa =
				Alignments.getPairwiseAlignment(query, target,
						PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
 
		System.out.println(psa);
                System.out.println(psa.getNumIdenticals());
                System.out.println(targetSeq.length());
                System.out.println(querySeq.length());
                System.out.println(psa.getLength());
                
                int st = 1;
                while (psa.hasGap(st))
                    st++;
                int end = psa.getLength();
                while (psa.hasGap(end))
                    end--;
                
                System.out.println(st);
                System.out.println(end);
                System.out.println(1.0 - ((double) psa.getNumIdenticals())/(end-st+1));
              
  
                
    }

}

