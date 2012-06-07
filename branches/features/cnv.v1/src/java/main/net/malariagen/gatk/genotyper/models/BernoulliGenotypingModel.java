package net.malariagen.gatk.genotyper.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import net.malariagen.gatk.genotyper.VariantPosteriors;
import net.malariagen.gatk.genotyper.GenotypingContext;
import net.malariagen.gatk.genotyper.models.GModel;
import net.malariagen.gatk.math.Beta;

import org.apache.commons.math.MathException;
import org.apache.commons.math.util.ContinuedFraction;
import org.apache.commons.math.util.FastMath;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.InferredGeneticContext;

@GModel(fullName = "HaploidClonalBernoulliModel", shortName = "HCBM", aliases = { "Bernoulli" }, description = "Model for genotyping clonal haploid samples")
public class BernoulliGenotypingModel extends AbstractGenotypingModel {

	public BernoulliGenotypingModel() {
	}

	public BernoulliGenotypingModel(GenotypingContext gc, double refProb) {
		this.varPriorProb = refProb;
		this.setGenotypingContext(gc);
		init();
	}

	@GParam(shortName = "vpp", fullName = "VarPriorProb", description = "Variant Prior Probability", offset = 0, required = true)
	protected double varPriorProb = 0.5;

	@GParam(shortName = "err", fullName = "ErrorRate", description = "Rate of error calls", offset = 1, required = false)
	protected double errorRate = 0;

	@GParam(fullName = "MaxHaplodiImbalance", shortName = "mhi", offset = 3, description = "Assuming these are a clonal haploid samples what is maximum fraction of reads that shall support the minor allele for the genotype to be considered good", required = false)
	public double maxHaploidImbalance = Double.NaN;

	private double refLogPrior;
	private double nefLogPrior;

	public void init() {
		super.init();
		if (Double.isNaN(varPriorProb))
			throw new IllegalArgumentException("the rate cannot be a NaN");
		if (varPriorProb < 0)
			throw new IllegalArgumentException("the rate cannot be negative");
		if (varPriorProb > 1)
			throw new IllegalArgumentException(
					"the rate cannot be greater than 1");
		refLogPrior = escale(1.0 - varPriorProb);
		nefLogPrior = escale(varPriorProb);
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof BernoulliGenotypingModel))
			return false;
		BernoulliGenotypingModel other = (BernoulliGenotypingModel) o;
		return Math.abs(other.varPriorProb - varPriorProb) < 0.0001;
	}

	@Override
	public int getGenotypeCount() {
		return 2;
	}

	@Override
	public List<Allele> getGenotypeAlleles(int genotype) {
		if (genotype == 0) {
			return Collections.singletonList(getGenotypingContext()
					.getReferenceAllele());
		} else if (genotype == 1) {
			return Collections.singletonList(getGenotypingContext().getAlleleCount() >= 2 ? getGenotypingContext()
					.getAllele(1) : Allele.create("N",false));
		} else {
			return Collections.singletonList(Allele.NO_CALL);
		}
	}

	@Override
	public double calculatePosteriors(AlignmentContext ac,
			VariantPosteriors dest, int sample) {

		ReadBackedPileup rb = ac.getBasePileup();

		if (this.getGenotypeCount() <= 1) {
			dest.setGenotypePosterior(sample, 0, 0);
			dest.setGenotypePosterior(sample, 1,
					(int) InferredGeneticContext.NO_NEG_LOG_10PERROR);
			dest.setGenotypeQuality(sample,
					InferredGeneticContext.NO_NEG_LOG_10PERROR);
			return 0;
		}
		List<Integer> refQuals = new ArrayList<Integer>();
		List<Integer> nefQuals = new ArrayList<Integer>();
		GenotypingContext gc = this.getGenotypingContext();
		byte refBase = gc.getReferenceAllele().getBases()[0];
		double refLogProb = refLogPrior;
		double nefLogProb = nefLogPrior;
		int nefCount = 0;
		int refCount = 0;
		double errorRateSum = 0;
		for (PileupElement e : rb) {
			if (e.isDeletion())
				continue;
			byte base = e.getBase();
			boolean isRef = base == refBase;
			byte qual = baseQual(e);
			if (isRef)
				refCount++;
			else if (!gc.hasAllele(base))
				continue;
			else
				nefCount++;
			double goodLogLk = NQUAL[qual];
			double errLogLk = qual;
			if (isRef) {
				refQuals.add((int)qual);
			}
			else {
				nefQuals.add((int)qual);
			}
			errorRateSum += unescale(errLogLk);
			double refLogLk = isRef ? goodLogLk : errLogLk;
			double nefLogLk = isRef ? errLogLk : goodLogLk;
			refLogProb += refLogLk;
			nefLogProb += nefLogLk;
		}
		double alpha = nefLogProb - refLogProb;
		double ratio = unescale(alpha);
		double ratio2 = unescale(-alpha);
		double errorRateAvg = errorRate > 0 ? errorRate : errorRateSum
				/ (refCount + nefCount);

		double postRefLogProb = escale(1.0 / (1.0 + ratio));
		double postNefLogProb = escale(1.0 / (1.0 + ratio2));
		if ((nefCount == 0 && refCount == 0)) {
			dest.setGenotypePosterior(sample, 0,
					(int) InferredGeneticContext.NO_NEG_LOG_10PERROR);
			dest.setGenotypePosterior(sample, 1,
					(int) InferredGeneticContext.NO_NEG_LOG_10PERROR);
			dest.setGenotypeQuality(sample, Double.NaN);
			dest.setGenotypeConfidence(sample, Double.NaN);
			dest.setAvgErrRate(sample, Double.NaN);
			return Double.NaN;
		}
		dest.setGenotypePosterior(sample, 0, postRefLogProb > 99999 ? 99999
				: postRefLogProb);
		dest.setGenotypePosterior(sample, 1, postNefLogProb > 99999 ? 99999
				: postNefLogProb);

		dest.setAvgErrRate(sample, escale(errorRateAvg));

		if (!Double.isNaN(maxHaploidImbalance)
				&& (postRefLogProb < postNefLogProb ? ((double) nefCount
						/ ((double) nefCount + refCount) > maxHaploidImbalance)
						: ((double) refCount / ((double) nefCount + refCount) > maxHaploidImbalance))) {
			dest.setGenotypeConfidence(sample, -99);
			return 0;
		}
		try {
			
			int x = postRefLogProb < postNefLogProb ? nefCount : refCount;

			double logPvalue;
			if (refCount > 0 && nefCount > 0) {
				logPvalue = Beta.log10(errorRateAvg, x + 1, refCount
						+ nefCount - x) * 10;
			} else if (postRefLogProb < postNefLogProb && nefCount == 0) {
				logPvalue = 0;
			} else if (postNefLogProb < postRefLogProb && refCount == 0) {
				logPvalue = 0;
			} else {
				// Should rarely happen.
				logPvalue = Beta.log10(errorRateAvg, x + 1.000001,
						refCount + nefCount - x + 0.000001) * 10;
			}
			dest.setGenotypeConfidence(sample, logPvalue < -99 ? -99
					: logPvalue);
			double qual = (postRefLogProb < postNefLogProb ? postNefLogProb
					- postRefLogProb : postRefLogProb - postNefLogProb);
			dest.setGenotypeQuality(sample, qual);
			return 1;
		} catch (Exception e1) {
			throw new GenotypingModelException(e1);
		}
	}

	private byte baseQual(PileupElement e) {
		int mq = e.getMappingQual();
		byte bq = e.getQual();
		byte result = mq < bq ? (byte) mq : bq;
		// 3 equal to 0.5 chances. 
		return result < 3 ? 3 : result;
	}

	public static double nqual(int q) {
		return -10 * Math.log10(1 - Math.pow(10, ((double) -q) / 10.0));
	}

	/**
	 * Precomputed Negative Q score table, indicating the Phred scaled
	 * probability that the base is a good call given its Q score (the error).
	 */
	protected static final double[] NQUAL;

	static {
		NQUAL = new double[100];
		for (int i = 0; i < NQUAL.length; i++) {
			NQUAL[i] = nqual(i);
		}
	}

	public double escale(double p) {
		return -10 * Math.log10(p);
	}

	public double unescale(double q) {
		return Math.pow(10, -q / 10.0);
	}

	@Override
	public GenotypePriors getGenotypePriors() {
		return new net.malariagen.gatk.genotyper.models.GenotypePriors(
				new double[] { varPriorProb, 1 - varPriorProb }, varPriorProb
						* (1 - varPriorProb) * 2);
	}

	public double calculateFrequentistVariantPhredQuality(
			Map<String, AlignmentContext> acs) {

		int nefCount = 0;
		int refCount = 0;
		double errorRateSum = 0;
		GenotypingContext gc = getGenotypingContext();
		byte refBase = gc.getReferenceAllele().getBases()[0];
		for (AlignmentContext ac : acs.values()) {
			ReadBackedPileup rb = ac.getBasePileup();
			for (PileupElement e : rb) {
				if (e.isDeletion())
					continue;
				byte b = e.getBase();
				boolean isRef = b == refBase;
				if (isRef)
					refCount++;
				else if (!gc.hasAllele(b))
					continue;
				else
					nefCount++;
				byte qual = baseQual(e);
				errorRateSum += unescale(qual);
			}
		}
		if (refCount + nefCount == 0)
			return Double.NaN;
		double errorRateAvg = errorRate > 0 ? errorRate : errorRateSum
				/ (refCount + nefCount);
		try {
			// int x = refCount < nefCount ? refCount : nefCount;
			// double logPvalue = log10RegularizedBeta(errorRate >= 0 ?
			// errorRate : errorRateAvg, x + 1, refCount + nefCount - x, 10e-15,
			// Integer.MAX_VALUE);

			int x = nefCount;
			double log10Pvalue;
			if (nefCount > 0 && refCount > 0) {
				log10Pvalue = Beta.log10(errorRateAvg, x + 1,
						refCount + nefCount - x);
			} else if (refCount == 0) {
				log10Pvalue = -999.99;
			} else {
				log10Pvalue = Beta.log10(errorRateAvg, x + 1.000001,
						refCount + nefCount - x + 0.000001);
			}

			double negLog10PError = -log10Pvalue;

			if (negLog10PError > 9999.99)
				negLog10PError = 9999.99;
			else if (negLog10PError < 0)
				negLog10PError = 0;
			return negLog10PError;
		} catch (RuntimeException e1) {
			throw new GenotypingModelException(e1);
		}
	}

	@Override
	public int getGenotypeIndex(List<Allele> alleles) {
		if (alleles.size() > 1 || alleles.size() < 1)
			return -1;
		Allele a = alleles.get(0);
		if (a == null || a.isNoCall() || a.isNull())
			return -1;
		byte b = a.getBases()[0];
		int i = getGenotypingContext().getAlleleIndex(b);
		if (i == 0)
			return 0;
		else if (i > 0)
			return 1;
		else
			return -1;
	}

}
