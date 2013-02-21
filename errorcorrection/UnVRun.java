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
public class UnVRun {
    public static void main(String[] args) throws IOException {



                        int mb = 1024*1024;
                        //Getting the runtime reference from system
                        Runtime runtime = Runtime.getRuntime();
                        System.out.println("##### Heap utilization statistics [MB] #####");
                        //Print used memory
                        System.out.println("Used Memory:"            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
                        //Print free memory
                        System.out.println("Free Memory:"            + runtime.freeMemory() / mb);
                        //Print total available memory
                        System.out.println("Total Memory:" + runtime.totalMemory() / mb);
                        //Print Maximum available memory
                        System.out.println("Max Memory:" + runtime.maxMemory() / mb);

                        int k = 25;
                        int f = 400;
                        int lt = 50;
                        int nIter = 3;
                        int mErPerc = 40;
                        int toClust = 1;

                        double minConcToPrint = 0.01;
                        int maxz = 1;
                        int kmin = k-15;
                        boolean toCalcHapl = true;
                        boolean toDelUncor = true;

//                        String dset_file = "Patient " + i + ", Time " + j + ".fas";
                        String dset_file = "unknownv.fasta";
                        String outdir = "results" + "(" + dset_file + ")";
                        DataSet ds = new DataSet(dset_file);
                        ds.setK(k);
                        ds.setFreqThr(f);
                        ds.setLenThr(lt);
                        ds.setMaxAllErrorsPerc(mErPerc);
                        ds.calculateKMersAndKCounts();
                        ds.PrintReadsDistrib(outdir);
/*                       Corrector cr = new Corrector(ds);
                        cr.setMaxz(maxz);
                        cr.setNIter(nIter);
                        cr.setToClust(toClust);
                        cr.setToRemoveAllUncorrect(toDelUncor);
                        cr.setKmin(kmin);
                        cr.setToFindHapl(toCalcHapl);
                        cr.run();
                        ds = cr.CorrectedReads();
                        ds.PrintCorrectedReads(dset_file+"_corrected.fas");
                        ds.PrintHaplotypes(dset_file+ "(k="+a[i][j]+")"+ "_haplotypes.fas");*/

		System.out.println("Finished!");

	}

}
