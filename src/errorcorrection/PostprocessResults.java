/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package errorcorrection;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author kki8
 */
public class PostprocessResults {
public static void main(String[] args) throws IOException, InterruptedException {


                int k = 25;
		int f = 70;
		int lt = 50;
		int nIter = 4;
		int mErPerc = 40;
		int toClust = 1;
                int toFindHapl = 0;
                int dominparamonenucl = 2;
                int nucldiffparam = 1;
                int dominparamgen = 2;
                int dominparampostpr = 2;
                boolean toPrintStat = true;
                int errorsseglen = (int) k/3;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = true;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???

//                String dset_file_name = "amp20.fas";
                String dset_file_name = "RL7_reversed.fas_corrected.fas";
//                String dset_file_name = "P7T6_reversed.fas_corrected.fas_haplotypes.fas.fas";
//                DataSet ds = new DataSet(dset_file_name);
                Corrector cr = new Corrector();
 //               ds.findHaplotypes();

//                ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                    cr.postprocessHaplotypes(dset_file_name,15,6.6,dominparampostpr);
                    System.out.println(dset_file_name);
                    cr.postprocessHaplotypesPairwise(dset_file_name +"_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
                    cr.postprocessHaplotypes(dset_file_name+"_postprocessed.fas_PostprocPair.fas",15,6.6,dominparampostpr);
                    cr.postprocessHaplotypesPairwise(dset_file_name + "_postprocessed.fas_PostprocPair.fas_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);

//		Corrector cr = new Corrector();
//                cr.printRevComp("Prep1-9.fas_corrected.fas_haplotypes.fas");
//                cr.postprocessHaplotypes("HOC_1_1b.fas_corrected.fas_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas_postprocessed.fas_PostprocPair.fas", 15, 6, 5);
//               cr.postprocessHaplotypesPairwise("HOC_1_1b.fas_corrected.fas_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair1.fas", 15, 6.6, 3, 10, nucldiffparam);


	}

}
