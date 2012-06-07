package net.malariagen.gatk.math;

import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.Arrays;


public class IntegerCounter {

	// values that get increase 
	private long count = 0;
	private long sum = 0;
	private long sqSum = 0;
	
	private int unknownCount = 0;

	// Before consolidating the counter
	// frequencies[i] contains the frequency of i. 
	
	private long[] frequencies = new long[1000]; // nice initial size.
	
	// short for values[value.length - 1];
	private int max = -1;
	// short for values[0];
	private int min = Integer.MAX_VALUE;
	// different depths encountered
	private int depthCount = 0;
	
	// 
	private ValueFrequencyArrayIntegerDistribution distribution;
	
	public IntegerCounter() {
		distribution = new ValueFrequencyArrayIntegerDistribution(0,0,0,new int[0],new long[0]);
	}
	
	@Deprecated
	public void applyIncrement(int depth) {
		addValue(depth);
	}
	
	public void addValue(int value) {
		
		if (value < 0) {
			unknownCount++;
			return;
		}
//		else if (depth == 0) {
//		    zeroDepthCount++;
//			return;
//		}
		distribution = null;
		count++;
		if (frequencies.length <= value) 
			frequencies = Arrays.copyOf(frequencies, value << 1);
		if (frequencies[value]++ == 0) depthCount++;
		sum += value;
		sqSum += value * value;
		if (value > max) max = value;
		if (value < min) min = value;
	}
	
	// Generate a distribution object from the current counter contents.
	ValueFrequencyArrayIntegerDistribution toDistribution() {
		if (distribution != null) 
			return distribution;
		
		int[] distValues = new int[depthCount];
		long[] distFrequencies = new long[depthCount];
		int nextIndex = 0;
		for (int v = min; v <= max; v++) {
			long f = frequencies[v];
			if (f == 0) continue;
			distValues[nextIndex] = v;
			distFrequencies[nextIndex++] = f;
		}
		return distribution = new ValueFrequencyArrayIntegerDistribution(count,sum,sqSum,distValues,distFrequencies);
	}

	// Add all the depth counts from the passed counter into this counter.
	public void applyCounter(IntegerCounter rhs) {
		if (frequencies.length <= rhs.max)
			frequencies = Arrays.copyOf(frequencies, rhs.max + 1);
		if (max < rhs.max) max = rhs.max;
		if (min > rhs.min) min = rhs.min;
		count += rhs.count;
		sum += rhs.sum;
		sqSum += rhs.sqSum;
		for (int i = rhs.min; i <= rhs.max; i++) {
			long rhsF = rhs.frequencies[i];
			if (rhsF == 0) continue;
			long f = this.frequencies[i];
			if (f == 0) depthCount++;
			frequencies[i] += rhsF;
		}
		distribution = null;
	}

	public void printReport(String indent, PrintWriter pw) {
		IntegerDistribution dist = toDistribution();
		NumberFormat NUMBER_FORMAT = NumberFormat.getInstance();
		NUMBER_FORMAT.setMaximumFractionDigits(4);
		pw.println(indent + count + " with known reference and non-zero coverage positions considered");
		if (unknownCount > 0) pw.println( indent + unknownCount + " unknown reference positions discarded");
		if (count == 0) return;
		pw.println(indent + sum + " total coverage (~ reads * read-length) ");
		pw.println(indent + NUMBER_FORMAT.format(dist.mean()) + " ± " + NUMBER_FORMAT.format(dist.standardDeviation()) + " reads mean coverage depth with cumulative probability " 
				+ NUMBER_FORMAT.format(dist.cumulativeProbability((int) Math.floor(dist.mean()))));
		pw.println(indent + NUMBER_FORMAT.format(dist.mode()) + " reads mode coverage depth with cumulative probability " 
				+ NUMBER_FORMAT.format(dist.cumulativeProbability(dist.mode())));
		pw.println(indent + 
				NUMBER_FORMAT.format(dist.minimum()) + ", " +
				NUMBER_FORMAT.format(dist.percentile(25)) + ", " +
				NUMBER_FORMAT.format(dist.median()) + ", " +
				NUMBER_FORMAT.format(dist.percentile(75)) + " and " +
				NUMBER_FORMAT.format(dist.maximum()) + 
				  " = min., Q1, median, Q3 and max. respectively");	
	}
	
		

}
