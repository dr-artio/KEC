/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package errorcorrection;

import errorcorrection.DataSet;
import java.io.IOException;

/**
 *
 * @author kki8
 */
public class RunEDAR {
        public static void main(String[] args) throws IOException {
	
//		String dset_file_name = "JFH1_titanium";
//		String dset_file_name = "test";
//		String dset_file_name = "test_small";
//		String dset_file_name = "test1";
//		String dset_file_name = "test2";
//		String dset_file_name = "Prep3-4";
//		String dset_file_name = "Prep3-1.fas";
//		String dset_file_name = "Prep3-2.fas";
//		String dset_file_name = "Prep3-3.fas";
//		String dset_file_name = "Prep2-1";
//		String dset_file_name = "Prep3-6";
//		String dset_file_name = "Prep3-11";
//		String dset_file_name = "Prep1-6.fas";
//		String dset_file_name = "Prep1-4.fas";
//		String dset_file_name = "Prep1-5.fas";
//		String dset_file_name = "Prep1-6";
//		String dset_file_name = "Prep1-7.fas";
//		String dset_file_name = "Prep1-8.fas";
//		String dset_file_name = "Prep1-9.fas";
//		String dset_file_name = "Prep1-10.fas";
		String dset_file_name = "Prep1-1.fas";
//		String dset_file_name = "Prep1-2.fas";
//		String dset_file_name = "Prep1-3.fas";
//		String dset_file_name = "Prep4-1";
//		String dset_file_name = "Prep4-2";
//		String dset_file_name = "Prep4-3";
//		String dset_file_name = "NS5A_JFA2001";
//		String dset_file_name = "HOC_27_1b";
//		String dset_file_name = "AHCV_01_1a";
//		String dset_file_name = "HCVNJ_7_1a";
//		String dset_file_name = "Kng003";
//		String dset_file_name = "Kng013";
//		String dset_file_name = "nhanes3hcvhvr1sequences";
//		String dset_file_name = "3.TCA.454Reads.fas";
//		String dset_file_name = "eps0.3_40reads.fas";
//		String dset_file_name = "R8_MID2_ID10";
//		String dset_file_name = "Patient 1, Time 1.fas";
//                String dset_file_name = "Patient 2, Time 6.fas";


		int k = 25;
		int f = 400;
		int lt = 50;
		int nIter = 3;
		int mErPerc = 40;
		int toClust = 1;

		double minConcToPrint = 0.2;
		int maxz = 1;
		int kmin = k-10;
		boolean toCalcHapl = true;

		String dset_file = dset_file_name;



		String outdir = "results" + "(" + dset_file + ")";
		DataSet ds = new DataSet(dset_file);
		ds.setK(k);
		ds.setFreqThr(f);
		ds.setLenThr(lt);
		ds.setMaxAllErrorsPerc(mErPerc);


//		ds.PrintUniqueReads(dset_file_name + "_unique.fas");

		Corrector cr = new Corrector(ds);
		cr.setMaxz(maxz);
		cr.setNIter(nIter);
		cr.setToClust(toClust);
		cr.setKmin(kmin);
		cr.setToFindHapl(toCalcHapl);
                cr.EDAR();
		ds = cr.CorrectedReads();
		ds.PrintCorrectedReads(dset_file_name+"_corrected_EDAR.fas");
		ds.PrintHaplotypes(dset_file_name+"_haplotypes_EDAR.fas");

		System.out.println("Finished!");

	}
}
