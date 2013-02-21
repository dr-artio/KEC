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
public class PyroNoiseRun {
     public static void main(String[] args) throws IOException, InterruptedException {

        int mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(rt.maxMemory()/mb);

		
//                int i = 2;
//                int j = 2;
//            String dset_file_name = "Even3.fas";
//              for (int i = 1; i <= 3; i++)
              {
                String dset_file_name = "Divergent.fas";
                int k = 25;
		int f = 75;
		int lt = 100;
		int nIter = 3;
		int mErPerc = 40;
		int toClust = 1;
                int toFindHapl = 0;
                int dominparamonenucl = 5;
                int nucldiffparam = 1;
                int dominparamgen = 25;
                int dominparampostpr = 25;
                int errorsseglen = 5;
                boolean toPrintStat = true;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = true;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???


		
/*                if (i <= 9)
                    dset_file_name = "Patient0" + i + ".fas";
                else
                    dset_file_name = "Patient" + i + ".fas";
//                String dset_file_name = "JFH1_titanium.fas";*/
		String dset_file = dset_file_name;

                File fl = new File(dset_file);
//                if (!fl.exists())
//                    continue;

		String outdir = "results" + "(" + dset_file + ")";
		DataSet ds = new DataSet(dset_file);
//                ds.fixDirection("ref.fas", 15, 6.6);
//                ds.PrintReads(dset_file + "_FixedDir.fas");


		ds.setK(k);
                ds.setFindErrorsSeglen(errorsseglen);
//		ds.setFreqThr(f);
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
//                ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                if (toFindHapl == 1)
                {
                    ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6, dominparampostpr);
//                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas",1,1);
                    cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
                }

		System.out.println("Finished!");
 //             }
                
               
                // SECOND RUN
                ds = null;
                cr = null;
                maxz = 0;
                nIter = 0;
                toPrintStat = true;
                toDelUncor = true;
                f = 100;
                toPostprocessHeur = false; //???
                dset_file_name = dset_file_name+"_corrected.fas";
		dset_file = dset_file_name;
                toFindHapl = 1;
                toCalcHapl = true;

		outdir = "results" + "(" + dset_file + ")";
		ds = new DataSet(dset_file);
		ds.setK(k);
                ds.setFindErrorsSeglen(errorsseglen);
//		ds.setFreqThr(f);
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
                ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                if (toFindHapl == 1)
                {
                    ds = null;
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6,dominparampostpr);
                    cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
                    System.out.println(dset_file_name);
                    cr.postprocessHaplotypesPairwise(dset_file_name +"_haplotypes.fas_postprocessed.fas_RevComp.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas",15,6.6,dominparampostpr);
                    cr.postprocessHaplotypesPairwise(dset_file_name + "_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
                }
//		cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");

		System.out.println("Finished!"); 

	}
     }
    
}
