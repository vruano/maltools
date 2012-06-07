package net.malariagen.gatk.math;

/**
 * Implementation of the Gamma functions used within net.malariagen
 * 
 * <p/>
 * Code is partly taken from Numeric Recipes 2.0 and adapted to Java.
 * 
 * @author Valentin Ruano-Rubio &lt;valentin.ruano@gmail.com&gt;
 */
public strictfp class Gamma {

	private static final double[] COF = { 76.18009172947146,
			-86.50532032941677, 24.01409824083091, -1.231739572450155,
			0.1208650973866179e-2, -0.5395239384953e-5 };

	private static final double LOG_10_INV = 1.0 / Math.log(10);
	
	public static double log(double xx) {
		int j;
		double x, y, tmp, ser;
		
		y = x = xx;
		tmp = x + 5.5;
		tmp -= (x + 0.5) * Math.log(tmp);
		ser = 1.000000000190015;
		for (j = 0; j < COF.length; j++)
			ser += COF[j] / ++y;
		return -tmp + Math.log(2.5066282746310005 * ser / x);
	}
	
	public static double log10(double xx) {
		return log(xx) * LOG_10_INV;
		
	}
}
