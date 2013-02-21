/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

import java.io.IOException;

/**
 *
 * @author kki8
 */
public class PostprocessShift {
    public static void main(String[] args) throws IOException, InterruptedException {


        int timepoints = 9;
        int p = 15;
//        int[] indexes = {3,4,5,6};
        int[] indexes = {1};
//        int i = 0;
        for (int i = 0; i < indexes.length; i++)
        {
 //               String dset_file_name = "P11T" + i + "_unique_meta_w4_no_ins_dels.fas";
                
                 int j = 0;
                if (indexes[i] <= 7)
                    j = 1;
                else
                    j = 2;
//		String dset_file_name = "S000" + (282 + i) + "_" + indexes[i] + ".R" + j + ".MID" + indexes[i] + " (single)_reversed.fas";
                
                String dset_file_name = "RL7_reversed.fas_corrected.fas_haplotypes.fas";
                
//                  String dset_file_name  = "S000" + (281+indexes[i]) + "_" + indexes[i] + "_unique_meta_w4_no_ins_dels.fas";
//                  String dset_file_name = "454ReadsMID" + indexes[i] + "_reversed.fas_corrected.fas_corrected.fas";
                Corrector cr = new Corrector();
//                cr.delGaps(dset_file_name);
                cr.postprocessHaplotypesNuclSwitch(dset_file_name,"KEC", 6.6, 6.6, 1);
//                cr.postprocessHaplotypesNuclSwitch(dset_file_name + "_delGaps.fas","KEC", 6.6, 6.6, 2);
//                cr.postprocessHaplotypesNuclSwitch(dset_file_name + "_delGaps.fas_postprocShift.fas","KEC", 6.6, 6.6, 2);
                
        }



	}

    
}
