/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package errorcorrection;

import java.io.File;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 *
 * @author kki8
 */
public class PostprocessShift {
    public static void main(String[] args) throws IOException, InterruptedException {


            String folder_name = "Postprocess1";
            File folder = new File(folder_name);
            System.out.println(folder.exists());
            File[] list_files = folder.listFiles();
            
        for (int i = 0; i < list_files.length; i++)
        {
                
                String dset_file_name = list_files[i].getPath();
                
//                String dset_file_name = "C13OK1_1_MID2_run1_reversed_unique.fas_PostprocShift.fas_Renamed.fas";
                
                Corrector cr = new Corrector();
//                cr.delGaps(dset_file_name);
//                cr.postprocessHaplotypesNuclSwitch(dset_file_name,"ET", 6.6, 6.6, 2);
                cr.postprocessHaplotypesNuclSwitch(dset_file_name,"ET", 6.6, 6.6, 2,5);
//                cr.postprocessAlignedHaplotypesNuclSwitch(dset_file_name,"ET", 2);

        }



	}

    
}
