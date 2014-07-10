/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

import ErrorCorrection.Corrector;
import java.io.File;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 *
 * @author kki8
 */
public class RenameReads {
     public static void main(String[] args) throws IOException, InterruptedException {


            String folder_name = "Postprocess";
            File folder = new File(folder_name);
            System.out.println(folder.exists());
            File[] list_files = folder.listFiles();
            
        for (int i = 0; i < list_files.length; i++)
        {
                
                String dset_file_name = list_files[i].getPath();
                
//                String dset_file_name = "C13OK1_1_MID2_run1_reversed_unique.fas_PostprocShift.fas_Renamed.fas";
                
                Corrector cr = new Corrector();

                String tag = "";
                StringTokenizer st = new StringTokenizer(dset_file_name,"_\\");
                st.nextToken();
                String s = st.nextToken();
                while (!s.equalsIgnoreCase("reversed"))
                {
                    tag += ("_" + s);
                    s = st.nextToken();
                }
                
                tag = tag.substring(1);
                cr.renameReads(dset_file_name, tag);
        }



	}
    
}
