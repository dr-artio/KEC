/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ErrorCorrection;

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
                int dominparamgen = 5;
                int dominparampostpr = 5;
                boolean toPrintStat = true;
                int errorsseglen = (int) k/3;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = true;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???
                
                
            String folder_name = "Postprocess";
            File folder = new File(folder_name);
            System.out.println(folder.exists());
            File[] list_files = folder.listFiles();
            
            double gapop = 15;
            double gapext = 6;
            
            for (int i = 0; i < list_files.length; i++)
            {

                    String dset_file_name = list_files[i].getPath();
                    
                    int ncorrections = 1;
                    while (ncorrections != 0)
                    {


    //                String dset_file_name = "P39_PT_sec2_unique_meta_w4_unique.fas_PostprocShift.fas";
    //                DataSet ds = new DataSet(dset_file_name);
                    Corrector cr = new Corrector();
    //                System.out.println(cr.findMinDist(dset_file_name1, dset_file_name2));
    //                ds.findHaplotypes();

    //                ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                    String idAlignAlg = "Muscle";
                    String idMethod = "KEC";
                        cr.postprocessHaplotypes(dset_file_name,gapop,gapext,dominparampostpr,idAlignAlg,idMethod);
                        System.out.println(dset_file_name);
                        ncorrections = cr.postprocessHaplotypesPairwise(dset_file_name +"_postprocessed.fas", gapop, gapext, dominparamonenucl, dominparamgen, nucldiffparam,idAlignAlg,idMethod);
                        dset_file_name = dset_file_name+"_postprocessed.fas_PostprocPair.fas";

                    }

            }
	}

}
