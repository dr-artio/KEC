package ErrorCorrection;

import java.io.File;
import java.io.IOException;

public class ErrorCorrection {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		
 /*               String folder_name = args[0];
                File folder = new File(folder_name);
                System.out.println(folder.exists());
                File[] list_files = folder.listFiles();*/
                
		String dset_file_name  = args[0];
		int k = Integer.parseInt(args[1]);
//		int lt = Integer.parseInt(args[2]);
//		int mErPerc = Integer.parseInt(args[2]);
		int nIter = Integer.parseInt(args[2]);
//		int toClust = Integer.parseInt(args[4]);
                int toFindHapl = Integer.parseInt(args[3]);
                int errorsseglen = Integer.parseInt(args[4]);
	

//                for (int i = 0; i < list_files.length; i++)
//              for (int i = 1; i <= 10; i++)
              {
                    
//                  String dset_file_name = list_files[i].getPath();
//                int k = 25;
		int f = 100;
		int lt = 50;
//		int nIter = 5;
		int mErPerc = 40;
		int toClust = 1;
//                int toFindHapl = 1;
                int dominparamonenucl = 5;
                int nucldiffparam = 1;
                int dominparamgen = 30;
                int dominparampostpr = 30; //30
                int minNReads = 900;
                boolean toPrintStat = true;
 //               int errorsseglen = (int) 25;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = true;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???


//		String dset_file_name = "Prep1-1.fas";
//                String dset_file_name = "JFH1_titanium.fas";
		String dset_file = dset_file_name;

                File fl = new File(dset_file);
//                if (!fl.exists())
//                    continue;

		String outdir = "results" + "(" + dset_file + ")";
		DataSet ds = new DataSet(dset_file);
		ds.setK(k);
//		ds.setFreqThr(f);
		ds.setLenThr(lt);
		ds.setMaxAllErrorsPerc(mErPerc);
                ds.setFindErrorsSeglen(errorsseglen);
		
		Corrector cr = new Corrector(ds);
//                cr.printRevComp("Prep1-9.fas_corrected.fas_haplotypes.fas");
//                cr.postprocessHaplotypesPairwise("HOC_34_1b.fas_corrected.fas_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas_postprocessed.fas_PostprocPair.fas", 15, 6.6, 3, 3, nucldiffparam);
		cr.setMaxz(maxz);
		cr.setNIter(nIter);
		cr.setToClust(toClust);
		cr.setKmin(kmin);
                cr.setToRemoveAllUncorrect(toDelUncor);
		cr.setToFindHapl(toCalcHapl);
                cr.setToPrintStat(toPrintStat);
                cr.setToPostprocessHeur(toPostprocessHeur);
                cr.setMinNReads(minNReads);
		cr.run();
		ds = cr.CorrectedReads();

/*                cr = new Corrector(ds);
		cr.setMaxz(maxz);
		cr.setNIter(1);
		cr.setToClust(toClust);
		cr.setKmin(kmin);
                cr.setToRemoveAllUncorrect(false);
		cr.setToFindHapl(toCalcHapl);
                cr.setToPrintStat(toPrintStat);
		cr.run();
		ds = cr.CorrectedReads();*/



		ds.PrintCorrectedReads(dset_file_name+"_corrected.fas");
/*                if (toFindHapl == 1)
                {
                    ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6, dominparampostpr);
//                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas",1,1);
                    cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
                }
				
		System.out.println("Finished!");*/


                // SECOND RUN
                maxz = 0;
                nIter = 0;
                toPrintStat = false;
                toDelUncor = true;
                toPostprocessHeur = false; //???
                dset_file_name = dset_file_name+"_corrected.fas";
		dset_file = dset_file_name;
//                toFindHapl = 1;

		outdir = "results" + "(" + dset_file + ")";
		ds = new DataSet(dset_file);
		ds.setK(k);
//		ds.setFreqThr(f);
		ds.setLenThr(lt);
		ds.setMaxAllErrorsPerc(mErPerc);
                ds.setFindErrorsSeglen(errorsseglen);
                minNReads = 1;


		cr = new Corrector(ds);
		cr.setMaxz(maxz);
		cr.setNIter(nIter);
		cr.setToClust(toClust);
		cr.setKmin(kmin);
                cr.setToRemoveAllUncorrect(toDelUncor);
		cr.setToFindHapl(toCalcHapl);
                cr.setToPrintStat(toPrintStat);
                cr.setToPostprocessHeur(toPostprocessHeur);
                cr.setMinNReads(minNReads);
		cr.run();
		ds = cr.CorrectedReads();

		ds.PrintCorrectedReads(dset_file_name+"_corrected.fas");
                ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
                if (toFindHapl == 1)
                {
/*                    dominparamgen = 5;
                    dominparampostpr = 5;
                    dominparamonenucl = 5;
                    nucldiffparam = 3;*/
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
