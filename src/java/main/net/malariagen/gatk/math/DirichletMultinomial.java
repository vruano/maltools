package net.malariagen.gatk.math;

import net.malariagen.gatk.math.Dirichlet;

public class DirichletMultinomial {
	
	private Dirichlet prior;
	
	public DirichletMultinomial(Dirichlet prior) {
		if (prior == null)
			throw new IllegalArgumentException("");
	}
	
	public double logPdf(int ... x) {
		if (x.length != prior.size())
			throw new IllegalArgumentException("x must have the same number of members as catetories in the prior");
		int N = 0;
		double A = 0;
		double result = 0;
		for (int i = 0; i < x.length; i++) {
			int n = x[i];
			double a = prior.alpha(i);
			N += n;
			A += a;
			result += Gamma.log(n + a) - Gamma.log(a);			
		}
		result += Gamma.log(A) - Gamma.log(N + A);
		return result;
	}
	

}
