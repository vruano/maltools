package net.malariagen.gatk.coverage;

import net.malariagen.gatk.math.IntegerDistribution;


public abstract class AbstractCoverageDistribution implements IntegerDistribution {

	@Override
	public double standardDeviation() {
		return Math.sqrt(variance());
	}

	@Override
	public double median() {
		return percentile(50);
	}

}
