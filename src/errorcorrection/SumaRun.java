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
            String folder_name = "HCV Outbreak batch 10-JF";
            File folder = new File(folder_name);
            System.out.println(folder.exists());
            File[] list_files = folder.listFiles();
        

//            int[] thresholds = {19,20,18,24,17,19,16,20,18,12,20,13,12,23,17,25,21,25,26,21,17,27,25,16,22,20,16,11,13,14,15,22,26,12,21,23,18,20,19,22,29,32,27,27,39,27,19,28,21,27,21,28};
//            int[] thresholds = {13,21,9,21,16,15,13,10,11,12,21,13,12,22,15,13,15,11,13,14,13,18,13,15,11,14,13,7,11,10,11,11,13,9,12,21,16,19,17,16,13,11,13,14,21,9,14,23,18,20,13,14};
//all_NH_110412            int[] thresholds = {18,20,21,13,16,19,16,24,14,20,24,27,28,23,16,9,13,14,19,20,33,12,18,26,12,23,12,24,21,26,26,21,17,13,19,16,22,7,18,12,12,11,17,14,15,13,21,26,12,13,16,18,15,13,14,27,14,22,25,27,27,21,27,1,19,29,22,17,21,28,18,13,21,19};
           int[] thresholds = {22,30,27,15,20,10,19,17,13,18,25};
            
            int nRun = 1;
            

            
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
