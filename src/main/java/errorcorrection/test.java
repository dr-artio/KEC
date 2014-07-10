/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ErrorCorrection;

import java.io.File;
import java.io.IOException;
import java.util.Stack;
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

        DataSet ds = new DataSet("prot.fas");
        ds.PrintUniqueReads("prot_unique.fas");
        Stack s = new Stack();
              
  
                
    }

}

