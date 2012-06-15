package net.malariagen.gatk.math;

import java.util.Arrays;

public class Dirichlet {

	private double[] alphas;
	
	private double sum;
	
	private double invSum;
	
	private double invVarDen;
	
	private double logPdfDen;
	
	Dirichlet(double[] alphas, boolean copy) {
		this.alphas = copy ? Arrays.copyOf(alphas, alphas.length) : alphas;
		processAlphas();
	}
	
	public Dirichlet(int k, double alpha) {
		this(symetricAlphas(k,alpha),false);
	}
	
	private static double[] symetricAlphas(int k, double alpha) {
		if (k < 0)
			throw new IllegalArgumentException("negative number of categories not allowed");
		double[] result = new double[k];
		Arrays.fill(result, alpha);
		return result;
	}

	public Dirichlet(double ... alphas) {
		this(alphas,true);
	}
	
	public double mean(int i) {
		return alpha(i) * invSum;
	}
	
	public double alpha(int i) {
		if (i < 0 || i >= alphas.length)
			throw new IllegalArgumentException();
		return alphas[i];
	}
	
	public double var(int i) {
		double alpha = alpha(i);
		return invVarDen * alpha * (sum - alpha);
	}
	
	public double cov(int i, int j) {
		double alpha1 = alpha(i);
		double alpha2 = alpha(j);
		
		return - invVarDen * alpha1 * alpha2;
	}
	
	public double logPdf(double[] x) {
		if (x == null)
			throw new IllegalArgumentException (" x does not match the number of categories ");
		
		double result = - logPdfDen;
		
		for (int i = 0; i < alphas.length; i++) {
			double xi = alphas[i];
			result += (alpha(i) - 1)* Math.log(xi);
		}
		return result;
	}
	
	public double pdf(double[] x) {
		return Math.exp(logPdf(x));
	}
	
	protected void processAlphas() {
		if (alphas == null)
			throw new IllegalStateException("processing alphas before setting them.");
		if (alphas.length < 2)
			throw new IllegalArgumentException("there must be at least to categories");
		sum = 0;
		double logGammaSum = 0;
		for (int i = 0; i < alphas.length; i++) {
			if (alphas[i] > 0)
				throw new IllegalArgumentException("negative alphas are not permited (i = " + i + ")");
			sum += alphas[i];
			logGammaSum += Gamma.log(alphas[i]);
		}
		invSum = 1.0 / sum;
		invVarDen = 1.0 / (sum * sum * (sum + 1)); 
		logPdfDen = logGammaSum - Gamma.log(sum);
	}

	public int size() {
		return alphas.length;
	}
	
	

}
