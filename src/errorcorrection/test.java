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
public class test {
    public static void main(String[] args) throws IOException, InterruptedException {

             int mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(rt.maxMemory()/mb);

		
//                int i = 9;
//                int j = 4;
            String folder_name = "Lilian";
            File folder = new File(folder_name);
            System.out.println(folder.exists());
            File[] list_files = folder.listFiles();
        

//            int[] thresholds = {19,20,18,24,17,19,16,20,18,12,20,13,12,23,17,25,21,25,26,21,17,27,25,16,22,20,16,11,13,14,15,22,26,12,21,23,18,20,19,22,29,32,27,27,39,27,19,28,21,27,21,28};
//            int[] thresholds = {13,21,9,21,16,15,13,10,11,12,21,13,12,22,15,13,15,11,13,14,13,18,13,15,11,14,13,7,11,10,11,11,13,9,12,21,16,19,17,16,13,11,13,14,21,9,14,23,18,20,13,14};
              int[] thresholds = {18,20,21,13,16,19,16,24,14,20,24,27,28,23,16,9,13,14,19,20,33,12,18,26,12,23,12,24,21,26,26,21,17,13,19,16,22,7,18,12,12,11,17,14,15,13,21,26,12,13,16,18,15,13,14,27,14,22,25,27,27,21,27,1,19,29,22,17,21,28,18,13,21,19};
//           int[] thresholds = {14,14,14,11};
            
            int nRun = 1;
            

              for (int i = 0; i < list_files.length; i++)
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
                
                 String dset_file_name = list_files[i].getPath();

		String dset_file = dset_file_name;

/*                File fl = new File(dset_file);
                if (!fl.exists())
                    continue;*/
                System.out.println(dset_file_name);

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
                    
                    ds.calculateKMersAndKCounts();
//                    ds.findFreqThresholdExponent();
                     ds.findFreqThresholdPoisson();
                    
                
	}
    }

}

