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
public class CalcRealErrorsStatRun {
    public static void main(String[] args) throws IOException, InterruptedException {

        
        int nSample = 8;
                DataSet ds1 = new DataSet("S" + nSample + ".fas");
//                ds1.findHaplotypes();
//                ds1.PrintHaplotypes("S1_hapl.fas");

//               String addr_pattern = "Unique_clones_by_timepoint.fas";
               String addr_pattern = "S" + nSample + "_strangeHapl.fas";
        
              
              DataSet ds = new DataSet(addr_pattern);
//              ds.PrintReadsSeparateFiles();
              int k = ds.reads.size();
              ds.PrintReadsSeparateFilesKeepName(addr_pattern);
              for (int j = 1; j <=1;j++)
              {
//                  String address = "RL" + j + ".fas";
                  String address = "S" + nSample + ".fas";
                  ds = new DataSet(address);
                  for (int i = 1; i <= k; i++)
                  {
    //                ds.PrintReadsSeparateFiles();
                    ds.calculateRealErrorsStat(addr_pattern + "_Read" + i + ".fas", 10, 15, 6.6);

                  }
              }
    }

}
