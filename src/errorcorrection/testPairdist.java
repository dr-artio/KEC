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
public class testPairdist {
    public static void main(String[] args) throws IOException, InterruptedException {

            DataSet ds = new DataSet("MSMhcvhvr1followup454sequences.fas");
            ds.printpairdistReads("pairdist.txt");
                
	}
    
}
