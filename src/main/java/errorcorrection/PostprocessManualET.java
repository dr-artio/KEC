/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author kki8
 */
public class PostprocessManualET {
    
    public static void main(String[] args) throws IOException, InterruptedException {


                String addr = "P15T2_unique_meta_w4_no_ins_dels.fas";
                Corrector cr = new Corrector();
                ArrayList<Integer> positions = new ArrayList<Integer>();
                ArrayList<Character> type = new ArrayList<Character>(); 
                ArrayList<Character> nucleotide = new ArrayList<Character>();
                
                positions.add(85);
                type.add('I');
                nucleotide.add('G');
                
 /*               positions.add(21);
                type.add('I');
                nucleotide.add('C');
                
                positions.add(85);
                type.add('I');
                nucleotide.add('G');
                
                positions.add(96);
                type.add('D');
                nucleotide.add('X');
                
                positions.add(134);
                type.add('I');
                nucleotide.add('G');*/
                
                
                cr.postprocessHaplotypesManualAllignedET(addr, positions, type, nucleotide);
                
	}
    
}
