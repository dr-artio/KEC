/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

import java.util.Comparator;

/**
 *
 * @author kki8
 */
public class ReadFreqComparator implements Comparator{

    public int compare(Object o1, Object o2) {
        Read r1 = (Read) o1;
        Read r2 = (Read) o2;
        return -r1.frequency + r2.frequency;
    }
    
}