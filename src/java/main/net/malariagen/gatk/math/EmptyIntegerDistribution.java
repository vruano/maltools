package net.malariagen.gatk.math;

public class EmptyIntegerDistribution implements IntegerDistribution {

	@Override
	public long count() {
		return 0;
	}

	@Override
	public double mean() {
		return Double.NaN;
	}

	@Override
	public double variance() {
		return Double.NaN;
	}

	@Override
	public double standardDeviation() {
		return Double.NaN;
	}

	@Override
	public double percentile(int pc) {
		return Double.NaN;
	}

	@Override
	public double quantile(double q) {
		return Double.NaN;
	}

	@Override
	public double median() {
		return Double.NaN;
	}

	@Override
	public int minimum() {
		return Integer.MIN_VALUE;
	}

	@Override
	public int maximum() {
		return Integer.MAX_VALUE;
	}

	@Override
	public double cumulativeProbability(int c) {
		return Double.NaN;
	}

	@Override
	public int mode() {
		return 0;
	}

}
