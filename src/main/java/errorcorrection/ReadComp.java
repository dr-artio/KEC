/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

/**
 *
 * @author kki8
 */
import java.util.Comparator;

public class ReadComp implements Comparator<Read>{
		public int compare(Read r1, Read r2) {
                    Integer d1 = r1.frequency;
                    Integer d2 = r2.frequency;
                    return -d1.compareTo(d2);
		}

}

