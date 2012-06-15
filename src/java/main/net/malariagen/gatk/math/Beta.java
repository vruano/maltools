package net.malariagen.gatk.math;

/**
 * Implementation of the Beta function used within net.malariagen 
 * 
 * <p/>
 * Notice tha 'strictfp' needed to be added in order to stop some inconsistent outcomes: despite being totally stateless in some (rare) instances i was returning 
 * different results for the same parameter set. So please do not remove. 
 * 
 * @author Valentin Ruano-Rubio &lt;valentin.ruano@gmail.com&gt;
 */
public strictfp class Beta {

	private static final double LOG_10_INV = 1.0 / Math.log(10);
	private static final int MAXIT = Integer.MAX_VALUE;
	private static final double EPS = 2.e-20;

	/**
	 * Returns the natural logarithm of the Beta function. 
	 * @param a alpha, must be greater than 0.
	 * @param b beta, must be greater than 0.
	 * @return return is not determined if {@code a <= 0 || b <= 0}; do not rely on it to be a {@link Double#NaN NaN}.
	 */
	public static double log(double a, double b) {
		return Gamma.log(a) + Gamma.log(b) - Gamma.log(a+b);
	}

	public static double log10(double a, double b) {
		return Gamma.log10(a) + Gamma.log10(b) - Gamma.log10(a+b);
	}
	
	public static double phred(double a, double b) {
		return -10 * log10(a,b);
	}
	
	// incomplete beta.
	public static double log(double x, double a, double b) {
		if (a <= 0.0 || b <= 0.0) return Double.NaN;
		if (x < 0.0 || x > 1.0) return Double.NaN;
		if (x == 0.0 || x == 1.0) return x;
		
		double bt = -log(a,b) + a * Math.log(x) + b * ((x < 0.001) ? Math.log1p(-x) : Math.log(1.0 - x));
		if (x < (a+1.0)/(a+b+2.0)) 
			return bt + logContinuedFraction(x,a,b) - Math.log(a);
		else {
			double y = Math.exp(bt + logContinuedFraction(1.0-x,b,a) - Math.log(b));
			return y < 0.001 ? Math.log1p(-y) : Math.log(1.0 - y);
		}
	}
	
	public static double log10(double x, double a, double b) {
		if (a <= 0.0 || b <= 0.0) return Double.NaN;
		if (x < 0.0 || x > 1.0) return Double.NaN;
		if (x == 0.0 || x == 1.0) return x;
		//if (a > SWITCH && b > SWITCH) return logByQuadrature(x,a,b);
		
		double bt = -log10(a,b) + a * Math.log10(x) + b * ((x < 0.001) ? Math.log1p(-x) * LOG_10_INV : Math.log10(1.0 - x));
		if (x < (a+1.0)/(a+b+2.0)) 
			return bt + log10ContinuedFraction(x,a,b) - Math.log10(a);
		else {
			double y = Math.pow(10,bt + log10ContinuedFraction(1.0-x,b,a) - Math.log10(b));
			return y < 0.001 ? Math.log1p(-y) * LOG_10_INV : Math.log10(1.0 - y);
		}
	}
	
	public static double phred(double x, double a, double b) {
		return -10.0D * Beta.log10(x,a,b);
	}

	private static double log10ContinuedFraction(double x, double a, double b) {
        return logContinuedFraction(x,a,b) * LOG_10_INV;
	}
	
	
	
	private static double logContinuedFraction(double x, double a, double b) {
		double qap,qam,qab,em,tem,d;
		double bz,bm=1.0,bp,bpp;
		double az=1.0,am=1.0,ap,app,aold = -1;
	

		qab=a+b;
		qap=a+1.0;
		qam=a-1.0;
		bz=1.0-qab*x/qap;
		for (int m=1;m<=MAXIT;m++) {
			em=(double) m;
			tem=em+em;
			d=em*(b-em)*x/((qam+tem)*(a+tem));
			ap=az+d*am;
			bp=bz+d*bm;
			d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
			app=ap+d*az;
			bpp=bp+d*bz;
			aold=az;
			am=ap/bpp;
			bm=bp/bpp;
			az=app/bpp;
			bz=1.0;
			if (Math.abs(az-aold) < (EPS*Math.abs(az))) return Math.log(az);
		}    	
		throw new RuntimeException(String.format("a or b too big, or MAXIT too small a=%g, b=%g, x=%g",a,b,x));
        
	}

}
