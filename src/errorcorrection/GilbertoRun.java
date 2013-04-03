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
public class GilbertoRun {
    public static void main(String[] args) throws IOException, InterruptedException {

             int mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(rt.maxMemory()/mb);

		
//                int i = 9;
//                int j = 4;
        
            int[] patients = {16};
            int[] indexes = {1,2,3,5,7,8,10,11,13,14};
            int[] thresholds = {25,15,14,16,18,13,28,12,24,17};
            

            
            int j=0;
//            for (int j = 0; j < patients.length; j++)
              for (int i = 0; i < indexes.length; i++)
              {

                int k = 25;
		int f = thresholds[i];
		int lt = 50;
		int nIter = 3;
		int mErPerc = 40;
		int toClust = 1;
                int toFindHapl = 0;
                int dominparamonenucl = 5;
                int nucldiffparam = 1;
                int dominparamgen = 25;
                int dominparampostpr = 25;
                int errorsseglen = 25;
                boolean toPrintStat = false;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = false;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???


		String dset_file_name = "MIC_2-" + indexes[i] + "_reversed.fas";
//                String dset_file_name = "P015W0_P.fas";
		String dset_file = dset_file_name;

                File fl = new File(dset_file);
                if (!fl.exists())
                    continue;

		String outdir = "results" + "(" + dset_file + ")";
/*		DataSet ds = new DataSet(dset_file);


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
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6, dominparampostpr);
//                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas",1,1);
                    cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
                }

                
		System.out.println("Finished!"); */


                // SECOND RUN
                DataSet ds = null;
                Corrector cr = null;
                maxz = 0;
                nIter = 0;
                toPrintStat = false;
                toDelUncor = true;
                f = thresholds[i];;
                toPostprocessHeur = false; //???
                dset_file_name = dset_file_name+"_corrected.fas";
		dset_file = dset_file_name;
                toFindHapl = 1;
                toCalcHapl = true;

		outdir = "results" + "(" + dset_file + ")";
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
                    ds = null;

                }
//		cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");

		System.out.println("Finished!"); 

            }
                
	}

}
