/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package errorcorrection;

import java.io.IOException;

/**
 *
 * @author kki8
 */
public class CalcRealErrorsStatRun {
    public static void main(String[] args) throws IOException, InterruptedException {


 //             DataSet ds1 = new DataSet("P7T6_reversed_junior.fas");
 //             ds1.printMostFreqReads("P7T6_reversed_junior.fas_20mostfreq.fas", 20);
               String address = "RL7_reversed.fas";
 //              addresses[1] = "VA16_1st_round_MID3_reversed.fas";
        
              String addr_pattern = "RL7_clones_rc.fas";
              DataSet ds = new DataSet(addr_pattern);
//              ds.PrintReadsSeparateFiles();
              int k = ds.reads.size();
              ds.PrintReadsSeparateFilesKeepName(addr_pattern);
//              for (int j = 0; j <2;j++)
              {
                  ds = new DataSet(address);
                  for (int i = 1; i <= k; i++)
                  {
    //                ds.PrintReadsSeparateFiles();
                    ds.calculateRealErrorsStat(addr_pattern + "_Read" + i + ".fas", 10, 15, 6.6);

                  }
              }
    }

}
