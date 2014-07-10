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
public class IBDMCRun {
    public static void main(String[] args) throws IOException, InterruptedException {

           int patient = 15;
           int[] indexes = {1,2,3,4,5,6,7,8,9};
           for (int i = 0; i < indexes.length; i++)
           {
               String dset_file_name = "P" + patient + "T" + indexes[i] + "_unique_meta_w4_no_ins_dels.fas_PostprocShift.fas_AA.fas";
               DataSet ds = new DataSet(dset_file_name, "KEC");
               ds.findHaplotypes();
               ds.PrintHaplotypes(dset_file_name + "_hapl.fas");
           }
              
/*           double s = 0;
           double q = 0;
           for (int i = 0; i < indexes.length; i++)
           {
               String addr = "P" + patient + "T" + indexes[i] + ".fas_corrected.fas_corrected.fas";
               DataSet ds = new DataSet(addr);
               int f1 = 0;
               int f2 = 0;
               for (Read r : ds.reads)
                   if (r.frequency == 1)
                       f1++;
                   else
                       f2++;
               double g1 = ((double) f1)/ds.reads.size();
               double g2 = ((double)f2)/ds.reads.size();
//               System.out.println(f1 + " " + f2);
               System.out.println(g1 + " " + g2);
                       
/*               File f = new File(addr);
               System.out.println(f.length());
               s+=f.length()/1024;    
               q+=ds.reads.size();*/
//           }
//           System.out.println(); 
//           System.out.println(s);
//           System.out.println(s/indexes.length);
//           System.out.println(q);
//           System.out.println(q/indexes.length);
                
	} 
                      
}
