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
public class SumaRun {
    public static void main(String[] args) throws IOException, InterruptedException {

             int mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(rt.maxMemory()/mb);

		
//                int i = 9;
//                int j = 4;
            String folder_name = "Pooling";
            File folder = new File(folder_name);
            System.out.println(folder.exists());
            File[] list_files = folder.listFiles();
        


           int[] thresholds = {3,5,6,8,4,6,4,7,7};
            
            int nRun = 2;
            

            
//            int j=0;
//            for (int j = 0; j < patients.length; j++)
              for (int i = 0; i < list_files.length; i++)
              {

                int k = 25;
		int f = thresholds[i];
		int lt = 50;
		int nIter = 3;
		int mErPerc = 50; // 40
		int toClust = 1;
                int dominparamonenucl = 5;
                int nucldiffparam = 1;
                int dominparamgen = 25;
                int dominparampostpr = 25;
                int errorsseglen = 25;
                boolean toPrintStat = false;

                int maxz = 3;
		int kmin = k-15;
                
                int toFindHapl = 0;
		boolean toCalcHapl = false;
                
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???
                

//		String dset_file_name = "S000" + (282 + i) + "_" + indexes[i] + ".R" + j + ".MID" + indexes[i] + " (single)_reversed.fas_correctedShiftReads.fas";
//                 String dset_file_name = "VA51_1st_round_MID4.fas";
//                 String dset_file_name = "454ReadsMID" + indexes[i] + "_reversed.fas";
                 String dset_file_name = list_files[i].getPath();
//                 String dset_file_name = "P" + indexes[i] + ".fas";

		String dset_file = dset_file_name;

/*                File fl = new File(dset_file);
                if (!fl.exists())
                    continue;*/

                if (nRun == 1)
                {
                    DataSet ds = new DataSet(dset_file);


                    ds.setK(k);
                    ds.setFindErrorsSeglen(errorsseglen);
                    ds.setFreqThr(f);
                    ds.setLenThr(lt);
                    ds.setMaxAllErrorsPerc(mErPerc);


                    Corrector cr = new Corrector(ds);


                    cr.setMaxz(maxz);
                    cr.setNIter(nIter);
                    cr.setToClust(toClust);
                    cr.setKmin(kmin);
                    cr.setToRemoveAllUncorrect(toDelUncor);
                    cr.setToFindHapl(toCalcHapl);
                    cr.setToPrintStat(toPrintStat);
                    cr.setToPostprocessHeur(toPostprocessHeur);
                    cr.run();
                    ds = cr.CorrectedReads();

                    ds.PrintCorrectedReads(dset_file_name+"_corrected.fas");
                    if (toFindHapl == 1)
                    {
                        ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
    //                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6, dominparampostpr);
    //                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas",1,1);
    //                    cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
                    }


                    System.out.println("Finished!"); 
                }


                if (nRun == 2)
                {
                    // SECOND RUN
                    DataSet ds = null;
                    Corrector cr = null;
                    maxz = 0;
                    nIter = 0;
                    toPrintStat = false;
                    toDelUncor = true;
                    f = thresholds[i];;
                    toPostprocessHeur = false; //???
//                    dset_file_name = dset_file_name+"_corrected.fas";
                    dset_file = dset_file_name;
                    toFindHapl = 1;
                    toCalcHapl = true;

                    ds = new DataSet(dset_file);
                    ds.setK(k);
                    ds.setFindErrorsSeglen(errorsseglen);
                    ds.setFreqThr(f);
                    ds.setLenThr(lt);
                    ds.setMaxAllErrorsPerc(mErPerc);


                    cr = new Corrector(ds);
                    cr.setMaxz(maxz);
                    cr.setNIter(nIter);
                    cr.setToClust(toClust);
                    cr.setKmin(kmin);
                    cr.setToRemoveAllUncorrect(toDelUncor);
                    cr.setToFindHapl(toCalcHapl);
                    cr.setToPrintStat(toPrintStat);
                    cr.setToPostprocessHeur(toPostprocessHeur);
                    cr.run();
                    ds = cr.CorrectedReads();

                    ds.PrintCorrectedReads(dset_file_name+"_corrected.fas");
                    if (toFindHapl == 1)
                    {
                        ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
//                                            cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6,dominparampostpr);
//                                            cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
//                                            System.out.println(dset_file_name);
//                                            cr.postprocessHaplotypesPairwise(dset_file_name +"_haplotypes.fas_postprocessed.fas_RevComp.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
//                                            cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas",15,6.6,dominparampostpr);
//                                            cr.postprocessHaplotypesPairwise(dset_file_name + "_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);

                    }
//                     cr.postprocessHaplotypesNuclSwitch(dset_file_name+"_haplotypes.fas","KEC", 6.6, 6.6, 2);
    //		cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");

                    System.out.println("Finished!"); 
                }

            }
                
	}

}
