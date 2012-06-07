package net.malariagen.gatk.math;

public interface IntegerDistribution {
	
	public long count();

	public double mean();
	
	public double variance();
	
	public double standardDeviation();
	
	public double percentile(int pc);
	
	public double quantile(double q);
	
	public double median();
	
	public int minimum();
	
	public int maximum();
	
	public double cumulativeProbability(int c);

	public int mode();
	
}
