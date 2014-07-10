/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ErrorCorrection;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author ava (artyomenkoav@gmail.com)
 */
public class Poisson {
    private static Double _lambda;
    private static Map<Double,Map<Integer,Double>> _cache = new HashMap<Double, Map<Integer, Double>>();
    
    /**
     * Calculate probability P(X<=k) == (lambda^k / k!)*e^(-lambda)
     * i. e. under assumption of Poisson distribution
     * @param k 
     *      integer number (k>=0)
     * @return
     *      probability P(X<=k) == (lambda^k / k!)*e^(-lambda)
     * @throws Exception 
     */
    public static double pPoiss(int k) throws Exception {
        return pPoiss(k, _lambda);
    }
    /**
     * Calculate probability P(X<=k) == (lambda^k / k!)*e^(-lambda)
     * i. e. under assumption of Poisson distribution
     * @param k 
     *      integer number (k>=0)
     * @param lambda
     *      distribution parameter (lambda == E(X) == Varience(X), lambda > 0)
     * @return
     *      probability P(X<=k) == (lambda^k / k!)*e^(-lambda)
     * @throws Exception 
     */
    public static double pPoiss(int k, Double lambda) throws Exception {
        if (k < 0 || lambda < 0) throw new Exception("Invalid arguments");
        if (k==0) return Math.exp(-lambda);
        if (isInCache(k, lambda)) return getFromCache(k, lambda);
        double result = pPois(k, lambda) + pPoiss(k-1,lambda);
        cache(k, lambda, result);
        return result;
    }
    
    /**
     * Calculate probability P(X==k) == (lambda^k / k!)*e^(-lambda)
     * i. e. under assumption of Poisson distribution
     * @param k 
     *      integer number (k>=0)
     * @param lambda
     *      distribution parameter (lambda == E(X) == Varience(X), lambda > 0)
     * @return
     *      probability P(X==k) == (lambda^k / k!)*e^(-lambda)
     * @throws Exception 
     */
    public static double pPois(int k, Double lambda) throws Exception {
        if (k <= 0 || lambda < 0) {
            throw new Exception("Parameter can't be less than zero");
        }
        Double currentValue = lambda;
        for (int i = 1; i<=k; i++){
            currentValue -= Math.log(lambda/i);
        }
        currentValue = Math.exp(-currentValue);
        
        return currentValue;
    }
    
    /**
     * Calculate probability P(X==k) == (lambda^k / k!)*e^(-lambda)
     * i. e. under assumption of Poisson distribution
     * @param k 
     *      integer number (k>=0)
     * @param stvar
     *      Map of values of stohastic variable
     * @return
     *      probability P(X==k) == (lambda^k / k!)*e^(-lambda)
     * @throws Exception 
     */
    public static double pPois(int k) throws Exception {
        return pPois(k, _lambda);
    }
    
    /**
     * Set _lambsa value.
     * @param stvar
     *      ordered array (index as a value of X)
     * @param start
     *      index of first value
     */
     public static void estimateLambdaPoisson(double[] stvar, int start){
         Map<Integer, Double> nstvar = new HashMap<Integer, Double>(stvar.length);
        Integer i = start;
        for (Double d : stvar){
            nstvar.put(i, d);
            i++;
        } 
        estimateLambdaPoisson(nstvar);
     }
     
     /**
     * Set _lambsa value.
     * @param stvar
     *      ordered array (index as a value of X)
     * @param start
     *      index of first value
     */
     public static void estimateLambdaPoisson(int[] stvar, int start){
        double[] nstvar = new double[stvar.length];
        for (int i = 0; i<stvar.length; i++) nstvar[i] = stvar[i];
        estimateLambdaPoisson(nstvar, start);
     }
     
     /**
     * Give estimation of lambda for discrete stohastic variable
     * and save to _lambda
     * @param stvar
     *      Map of values of stohastic variable
     * @return 
     *      lambda (Parameter for Poisson distribution)
     */
    public static void estimateLambdaPoisson(Map<Integer, Double> stvar){
        double sum = 0;
        int size = 0;
        for (int i : stvar.keySet()){
            sum += stvar.get(i) * i;
            size += stvar.get(i);
        }
        _lambda = sum / size;
        System.out.println(_lambda);
    }
    
    private static boolean isInCache(int k, double lambda){
        return _cache.containsKey(lambda) && ((Map<Integer, Double>)_cache.get(lambda)).containsKey(k);
    }
    
    private static Double getFromCache(int k, double lambda){
        return ((Map<Integer, Double>)_cache.get(lambda)).get(k);
    }
    
    private static void cache(int k, double lambda, double value){
        Map<Integer, Double> mapForLambda;
        if (!_cache.containsKey(lambda)){
            mapForLambda = new HashMap<Integer, Double>();
            _cache.put(lambda, mapForLambda);
        } else mapForLambda = _cache.get(lambda);
        mapForLambda.put(k, value);
    }
}
