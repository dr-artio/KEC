/*
 * Copyrights
 */
package ErrorCorrection;

/**
 * This class is for making dynamic outpus
 * for Console applications. Works with 
 * System.out stream. Usual ussage for making
 * progress bar or indicateng steps passed
 * 
 * @author alexander(artyomenkoav@gmail.com)
 */
public class DynamicOut {
    private static String _previous = "";
    
    /**
     * This method is print string into 
     * Sustem.out and hide it before next print
     * @param s 
     *      String temporary displayed on console
     */
    public static void printStep(String s){
        hideCurrent();
        _previous = s;
        System.out.print(s);
    }
    
    /**
     * Call this method after cycle where used 
     * PrintStep() method
     */
    public static void finishSteps(){
        _previous = "";
        System.out.println();
    }
    private static void hideCurrent(){
        StringBuilder sb = new StringBuilder();
        int length = _previous.length();
        for (int i = 0; i <= length; i++) sb.append("\b");
        System.out.print(sb.toString());
    }
}
