/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

import java.io.File;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 *
 * @author kki8
 */
public class CalcStat {
    public static void main(String[] args) throws IOException, InterruptedException {

        int mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(rt.maxMemory()/mb);
        
        String folder_name = "Palmer";
        File folder = new File(folder_name);
        System.out.println(folder.exists());
        File[] list_files = folder.listFiles();
        
        int rename = 0;
        
        if (rename == 1)
        {
            for (int i = 0; i < list_files.length; i++)
            {
                 StringTokenizer st = new StringTokenizer(list_files[i].getPath(),"\\");
                 String s0 = st.nextToken();
                 String s1 = st.nextToken();
                 String s2 = "";
                 if (i < 10)
                     s2 = "0" + i;
                 else
                     s2 = "" + i;
                 File f = new File(s0 + "\\" + s2 + "_" + s1);
                 list_files[i].renameTo(f);
            }
            list_files = folder.listFiles();
        }

		int[] indexes = {2,3,4,5,7,8,9,10};

              for (int i = 0; i < list_files.length; i++)
              {

                int k = 25;
		int f = 63;
		int lt = 50;
		int nIter = 3;
		int mErPerc = 40;
		int toClust = 1;
                int toFindHapl = 0;
                int dominparamonenucl = 5;
                int nucldiffparam = 1;
                int dominparamgen = 25;
                int dominparampostpr = 25;
                int errorsseglen = 20;
                boolean toPrintStat = true;

                int maxz = 3;
		int kmin = k-15;
		boolean toCalcHapl = false;
                boolean toDelUncor = false;
                boolean toPostprocessHeur = false; //???


                 String dset_file_name = list_files[i].getPath();
 
		String dset_file = dset_file_name;

                File fl = new File(dset_file);
//                if (!fl.exists())
//                    continue;

                StringTokenizer st = new StringTokenizer(dset_file,"\\");
                st.nextToken();
                String s1 = st.nextToken();
		String outdir = "results" + "(" + s1 + ")";
//                String outdir = "results" + "(" +  dset_file + ")";
		DataSet ds = new DataSet(dset_file);
//                ds.fixDirection("ref.fas", 15, 6.6);
//                ds.PrintReads(dset_file + "_FixedDir.fas");


		ds.setK(k);
                ds.setFindErrorsSeglen(errorsseglen);
		ds.setFreqThr(f);
		ds.setLenThr(lt);
		ds.setMaxAllErrorsPerc(mErPerc);
                
                ds.calculateKMersAndKCounts();
                ds.PrintReadsStat(outdir);
		
	}

    

    }
}
