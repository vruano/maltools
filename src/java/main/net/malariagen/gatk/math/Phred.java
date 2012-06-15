package net.malariagen.gatk.math;

public strictfp class Phred {

	public static double doubleValue(double p) {
		return -10.0D * Math.log10(p);
	}
	
	public static String stringValue(double p, int decimals) {
		if (decimals < 0)
			return "" + p;
		if (decimals == 0 ) return String.format("%d",doubleValue(p));
		return String.format("%" + decimals + "f",doubleValue(p));
	}
	
	public static String stringValue(double p) {
		return stringValue(p,2);
	}
	
	public static double prob(double d) {
		return Math.pow(10, - d * 0.1D);
	}
}
