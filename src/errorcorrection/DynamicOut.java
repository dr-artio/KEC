/*
 * Copyrights
 */
package errorcorrection;

/**
 * This class is for making dynamic output
 * for Console applications. Works with 
 * System.out stream. Usual usage for making
 * progress bar or indicating steps passed
 * 
 * @author alexander(artyomenkoav@gmail.com)
 */
public class DynamicOut {
    private static String _previous = "";
    
    /**
     * This method is print string into 
     * System.out and hide it before next print
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
