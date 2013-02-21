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
public class KNGRun {
    public static void main(String[] args) throws IOException, InterruptedException {

        //                int i = 5;
              for (int i = 1; i <= 9; i++)
              {

                int k = 25;
		int f = 318;
		int lt = 50;
		int nIter = 4;
		int mErPerc = 40;
		int toClust = 1;
                int toFindHapl = 1;
                int dominparamonenucl = 4;
                int nucldiffparam = 1;
                int dominparamgen = 30;
                int dominparampostpr = 40;
                boolean toPrintStat = false;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = true;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = true; //???


		String dset_file_name = "Kng00" + i + ".fas";
		String dset_file = dset_file_name;

		String outdir = "results" + "(" + dset_file + ")";
		DataSet ds = new DataSet(dset_file);
		ds.setK(k);
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
                if (toFindHapl == 1)
                {
                    ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6,dominparampostpr);
//                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas",1,1);
                }
		cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
		System.out.println("Finished!");


                // SECOND RUN
                maxz = 0;
                nIter = 0;
                toPrintStat = true;
                toDelUncor = true;
                toPostprocessHeur = false; //???
                dset_file_name = dset_file_name+"_corrected.fas";
		dset_file = dset_file_name;

		outdir = "results" + "(" + dset_file + ")";
		ds = new DataSet(dset_file);
		ds.setK(k);
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
                if (toFindHapl == 1)
                {
                    ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6,dominparampostpr);
                    cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
                    System.out.println(dset_file_name);
                    cr.postprocessHaplotypesPairwise(dset_file_name +"_haplotypes.fas_postprocessed.fas_RevComp.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
                }
//		cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");

		System.out.println("Finished!");

            }
        

	}


}
