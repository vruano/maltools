package net.malariagen.gatk.genotyper.models;


public class GenotypePriors implements org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors {

	protected double heterozygosity;
	
	protected double[] priors;
	
	GenotypePriors(double[] priors, double heterozygosity) {
		this.heterozygosity = heterozygosity;
		this.priors = priors;
	}
	
	@Override
	public double[] getPriors() {
		return priors;
	}

	@Override
	public double getHeterozygosity() {
		return heterozygosity;
	}

	@Override
	public boolean validate(boolean throwException) {
		return true;
	}

}
