package net.malariagen.gatk.genotyper.models;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import net.malariagen.gatk.genotyper.VariantPosteriors;
import net.malariagen.gatk.genotyper.GenotypingContext;
import net.malariagen.utils.NucleotideIUPAC;

import org.apache.commons.math.util.MathUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import scala.actors.threadpool.Arrays;

@GModel(fullName = "DiscretePolyploidMixtureModel", shortName = "DPMM", aliases = { "DiscreteMixture" }, description = "Model for genotyping clonal haploid samples")
public class DiscreteMixtureGenotypingModel extends AbstractGenotypingModel {

	@GParam(shortName = "prior", fullName = "PriorDensity", description = "file containing the prior probability densities", offset = 0, required = true)
	protected File prior = null;

	@GParam(shortName = "err", fullName = "ErrorRate", description = "Rate of error calls", offset = 1, required = false)
	protected double errorRate = 0;
	
	private int genotypeCount = 0;
	
	private double[] priorProbs;
	
	private double[] oneMinusPriorProbs;
	
	private double[] priorFreqs;
	
	private GenotypePriors genotypePriors;
	
	public DiscreteMixtureGenotypingModel() {
		prior = null;
		errorRate = 0;
	}
	
	public DiscreteMixtureGenotypingModel(File prior, double errorRate) {
		if (Double.isNaN(errorRate))
			throw new IllegalArgumentException("invalid error rate, must not be NaN");
		if (errorRate < 0)
			throw new IllegalArgumentException("invalid error rate, must not be negative");
		if (Double.isInfinite(errorRate))
			throw new IllegalArgumentException("invalid error rate, must not be infinite");
		if (prior == null)
			throw new NullPointerException();
		this.prior = prior;
		this.errorRate = errorRate;
	}
	
	@Override
	void init() {
		super.init();
		loadPrior();
		buildGenotypePriors();
		GenotypingContext gc = this.getGenotypingContext();
		int alleleCount = gc.getAlleleCount();
		
		if (genotypeCount != (1 << alleleCount) - 1) {
			throw new GenotypingModelException(String.format("mismatch between prior genotypes (%d) and number of alleles (%d)",genotypeCount,(1 << alleleCount) - 1));
		}
	}
	
	private void loadPrior() {
		if (priorProbs != null)
			return;
		loadPriorFile();
	}
	
	private void loadPriorFile() {
		if (prior == null)
			throw new IllegalStateException("prior model parameter not set before executing init");
		if (!prior.exists())
			throw new GenotypingModelException(String.format("the prior file '%s' is not reachable or does not exists",prior.toString()));
		if (!prior.isFile())
			throw new GenotypingModelException(String.format("the prior file '%s' is in fact not a regular file",prior.toString()));
		if (!prior.canRead())
			throw new GenotypingModelException(String.format("the prior file '%s' cannot be read",prior.toString()));
		BufferedReader r;
		
		try {
			r = new BufferedReader(new FileReader(prior));
			PriorHeader ph = new PriorHeader(r.readLine());
			if (ph.alleleCount > 2)
				throw new GenotypingModelException("more than two alleles not yet supported");
			
			genotypeCount = (1 << ph.alleleCount) - 1;
			
			String line;
			priorProbs = new double[ph.categoryCount];
			priorFreqs = new double[ph.categoryCount];
			oneMinusPriorProbs = new double[ph.categoryCount];
			Integer[] indexes = new Integer[ph.categoryCount];
			double probTotal = 0;
			int idx = 0;
			while ((line = r.readLine()) != null) {
				indexes[idx] = idx;
				probTotal += ph.readLine(line,priorProbs,priorFreqs,oneMinusPriorProbs,idx);
				double prob = priorProbs[idx];
				if (prob < 0)
					throw new GenotypingModelException("illeal prior probability");
				idx++;
			}
			double missingProb = 1 - probTotal;
			if (Math.abs(missingProb) > 0.00001 * ph.categoryCount)
				throw new GenotypingModelException("prior probability leakage (1 != " + missingProb + ")");
			Arrays.sort(indexes, new Comparator<Integer>() {

				@Override
				public int compare(Integer arg0, Integer arg1) {
					if (priorFreqs[arg0] < priorFreqs[arg1])
						return -1;
					else if (priorFreqs[arg1] < priorFreqs[arg0])
						return 1;
					else 
						return arg0 < arg1 ? -1 : 1;
				}});
			double[] newPriorFreqs = new double[ph.categoryCount];
			double[] newPriorProbs = new double[ph.categoryCount];
			for (int i = 0; i < ph.categoryCount; i++) {
				newPriorFreqs[i] = priorFreqs[indexes[i]];
			    newPriorProbs[i] = priorProbs[indexes[i]];
			    oneMinusPriorProbs[i] = oneMinusPriorProbs[indexes[i]];
			}
			priorFreqs = newPriorFreqs;
			priorProbs = newPriorProbs;
			if (priorFreqs[0] != 0) {
				Logger.getLogger(this.getClass()).warn("there is no category with frequency 0 in prior");
			}
			if (priorFreqs[priorFreqs.length - 1] != 1) {
				Logger.getLogger(this.getClass()).warn("there is no category with frequency 1 in prior");
			}
			
		} catch (IOException e) {
			throw new GenotypingModelException(String.format("unexpected I/O exception encountered while loading prior file '%s':",prior.toString()),e);
		}
	}

	@Override
	public int getGenotypeCount() {
		return genotypeCount;
	}

	@Override
	public GenotypePriors getGenotypePriors() {
		if (genotypePriors == null)
			throw new IllegalStateException("prior requested before set-up");
		return genotypePriors;
	}
	
	public void buildGenotypePriors() {
		double[] priors = new double[this.getGenotypeCount()];
		if (priors.length == 1)
			priors[0] = 0;
		else if (priors.length == 3) {
			priors[0] =  (priorFreqs[0] <= 0.0000001) ? priorProbs[0] : Double.POSITIVE_INFINITY; 
			priors[1] =  (priorFreqs[priorFreqs.length - 1] >= 0.999999) ? priorProbs[priorProbs.length - 1] : Double.POSITIVE_INFINITY;
			priors[2] = (Double.isInfinite(priors[0]) && Double.isInfinite(priors[1])) ? 0 : 
				C1P * Math.log1p(- (Double.isInfinite(priors[0]) ? 0 : Math.pow(10, - priors[0] *0.1)) - 
					(Double.isInfinite(priors[1]) ? 0 : Math.pow(10, - priors[1] * 0.1)));
		}
		else
			throw new GenotypingModelException("allele count > 2 unsupported, but genotype count is neither 1 or 3");

		double hetero = 0;
		for (int i = 0; i < priorFreqs.length;i++) {
			double p = Math.pow(10, - 0.1 * priorProbs[i]);
			double h = 2 * priorFreqs[i] * ( 1 - priorFreqs[i]);
		    hetero += p * h;
		}
		
		for (int i = 0; i < priors.length; i++)
			priors[i] = Math.pow(10, - 0.1 * priors[i]);
		
		this.genotypePriors = new net.malariagen.gatk.genotyper.models.GenotypePriors(priors,hetero);
	}

	@Override
	public List<Allele> getGenotypeAlleles(int genotype) {
		if (genotype > getGenotypeCount())
			throw new IllegalArgumentException();
		if (genotype < 0 || genotype >= this.getGenotypeCount()) { 
			return Collections.singletonList(Allele.NO_CALL);
		}
		
		GenotypingContext gc = getGenotypingContext();
		List<Allele> allAlleles = gc.getAlleleList();
		int alleleCount = allAlleles.size();
		genotype++;
		int mask = 1;
		int resultAlleleCount = 0;
		int firstFound = -1;
		for (int i = 0; i < alleleCount; i++) {
			if ((genotype & mask) != 0) {
				if (firstFound < 0)
					firstFound = i;
				resultAlleleCount++;
			}
			mask <<= 1;
		}
		if (resultAlleleCount == 0)
			throw new GenotypingModelException("Number of alleles cannot be 0");
		// Quick solution for a single allele genotypes.
		if (resultAlleleCount == 1)
			return Collections.singletonList(gc.getAllele(firstFound));
		int idx = 0;
		List<Allele> result = new ArrayList<Allele>(resultAlleleCount);
		while (genotype != 0) {
			if ((genotype & 1) != 0)
				result.add(allAlleles.get(idx));
			genotype >>= 1;
			idx++;
		}
		return result;
	}

	@Override
	public int getGenotypeIndex(List<Allele> alleles) {
		if (alleles.size() != 1)
			return -1;
		Allele a = alleles.get(0);
		byte b = a.getBases()[0];
		GenotypingContext gc = getGenotypingContext();
		int aIdx = gc.getAlleleIndex(b);
		if (aIdx != -1)
			return aIdx;
		NucleotideIUPAC n = NucleotideIUPAC.fromBases(gc.getReferenceAllele().getBases()[0],b);
		if (n.byteValue() == b)
			return gc.getAlleleCount();
		else
			return -1;
	}

	@Override
	protected double calculatePosteriors(AlignmentContext ac,
			VariantPosteriors dest, int sample) {
	
		int refCount = 0;
		int nefCount = 0;
		double errorSum = 0;
		GenotypingContext gc = getGenotypingContext();
		Allele refAllele = gc.getReferenceAllele();
		byte refByte = refAllele.getBases()[0];
		for (PileupElement pe : ac.getBasePileup()) {
			byte b = pe.getBase();
			byte bq = pe.getQual();
			int mq = pe.getMappingQual();
			int q = (bq < mq) ? bq : mq;
			q = (q < 3) ? 3 : q; // 3 is 0.5 chance, the absolute minimum.
			if (b == refByte) refCount++; else nefCount++;
			errorSum += Math.pow(10,- 0.1  * q);
		}
		int totalCount = refCount + nefCount;
		double errorAvgRate = errorRate > 0 ? errorRate : errorSum / totalCount; 
		double[] log10Numerators =  new double[priorProbs.length];
		double pbc = C1P * MathUtils.binomialCoefficientLog(totalCount,nefCount);
		double minLogExp = Double.POSITIVE_INFINITY;
		double maxLogExp = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < priorProbs.length; i++) { 
		    double pNef = priorFreqs[i] + errorAvgRate - 2 * priorFreqs[i] * errorAvgRate;
		    double phredNef = phredOf(pNef);
		    double phredRef = phredOf1p(pNef);
		    log10Numerators[i] = - 0.1 * (priorProbs[i] + pbc + refCount * phredRef + nefCount * phredNef);
		    if (log10Numerators[i] < minLogExp)
		    	minLogExp = log10Numerators[i];
		    if (log10Numerators[i] > maxLogExp) {
		    	maxLogExp = log10Numerators[i];
		    }
		}
		double avgLogExp = (maxLogExp + minLogExp) / 2.0;
		double c = Math.max(maxLogExp - 300,avgLogExp);
		double gt0SumMC = 0, gt1SumMC = 0, gt2SumMC = 0;
		double allSumMC = 0, sumFreq = 0, sqSumFreq = 0;
		int gt1NumeratorCount = 0, gt2NumeratorCount = 0, gt0NumeratorCount = 0;
		double gt1Log10Sum = 0, gt2Log10Sum = 0, gt0Log10Sum = 0;
		double[] linearNumeratorMinusC = log10Numerators; // alias, will overwrite.
		for (int i = 0; i < priorProbs.length; i++) {
			if (priorFreqs[i] < 0.000001) 
				gt0Log10Sum += log10Numerators[i];
			else if (priorFreqs[i] > 0.999999) 
				gt1Log10Sum += log10Numerators[i];
			else
				gt2Log10Sum += log10Numerators[i];
			linearNumeratorMinusC[i] = Math.pow(10,log10Numerators[i] - c);			
			double freq =  priorFreqs[i] * linearNumeratorMinusC[i];
			sumFreq += freq;
			sqSumFreq += priorFreqs[i] * freq; 
			allSumMC += linearNumeratorMinusC[i];
			if (priorFreqs[i] < 0.000001) {
			    gt0SumMC += linearNumeratorMinusC[i];
			    gt0NumeratorCount++;
			}
			else if (priorFreqs[i] > 0.999999) {
			    gt1SumMC += linearNumeratorMinusC[i];
			    gt1NumeratorCount++;
			}
			else { 
			    gt2SumMC += linearNumeratorMinusC[i];
			    gt2NumeratorCount++;
			}
		}
		double avgFreq = sumFreq / allSumMC;
		double sdFreq = Math.sqrt(Math.abs(avgFreq * avgFreq - sqSumFreq / allSumMC));
	    double gt0PrePostMC,gt1PrePostMC,gt2PrePostMC;
		double allPhredSum = - 10 * (Math.log10(allSumMC) + c);
		gt0PrePostMC = gt0NumeratorCount == 1 ? - 10 * (gt0Log10Sum) - allPhredSum  : - 10 * (Math.log10(gt0SumMC) + c) - allPhredSum;
		gt1PrePostMC = gt1NumeratorCount == 1 ? - 10 * (gt1Log10Sum) - allPhredSum  : - 10 * (Math.log10(gt1SumMC) + c) - allPhredSum;
		gt2PrePostMC = gt2NumeratorCount == 1 ? - 10 * (gt2Log10Sum) - allPhredSum  : - 10 * (Math.log10(gt2SumMC) + c) - allPhredSum;
	
		dest.setSamplePosteriors(sample, gt0PrePostMC, gt1PrePostMC, gt2PrePostMC);
		int genotype = setGenotypeQuality(dest, sample, gt0PrePostMC, gt1PrePostMC, gt2PrePostMC);
	    dest.setAvgErrRate(sample, C1P * Math.log(errorAvgRate));
		dest.setAttribute(sample,"PF",avgFreq);
		dest.setAttribute(sample,"SF",sdFreq);
		setBayesianPriors(dest, sample, genotype, gt0PrePostMC, gt1PrePostMC, gt2PrePostMC);
		return 1;
	}
	
	private void setBayesianPriors(VariantPosteriors dest, int sample, int genotype, double ... posteriors) {
		GenotypePriors gp = getGenotypePriors();
		double[] linearPriors = gp.getPriors();
		int count = linearPriors.length;
	    double[] priors = new double[count];
	    for (int i = 0; i < count; i++)
	    	priors[i] = -10 * Math.log10(linearPriors[i]);
	    int outputCount =(count * (count + 1)) >> 1;
	    if (outputCount == 0)
	    	return;
	    StringBuffer sb = new StringBuffer(outputCount * 10);
	    double worstBfactor = Double.POSITIVE_INFINITY;
	    for (int i = 0; i < count; i++) 
	    	for (int j = i + 1; j < count; j++) {
	    		double bf = - (posteriors[i] + priors[j] - posteriors[j] - priors[i]);
	    		if (i == genotype) {
	    			if (bf < worstBfactor)
	    				worstBfactor = bf;
	    		}
	    		else if (j == genotype) {
	    			if (- bf < worstBfactor)
	    				worstBfactor = - bf;
	    		}
	    		sb.append(String.format("%.2f", bf)).append(',');
	    	}
	    sb.setLength(sb.length() - 1);
//	    System.err.println(String.format("%.2f %.2f %.2f",priors[0],priors[1],priors[2]));
	    dest.setAttribute(sample, "BF", sb.toString());
	    if (worstBfactor != Double.POSITIVE_INFINITY) dest.setGenotypeConfidence(sample, worstBfactor);
	}

	private int setGenotypeQuality(VariantPosteriors dest, int sample,
			double gt0Post, double gt1Post, double gt2Post) {

		if (gt2Post <  gt1Post) {
			if (gt2Post < gt0Post) {
				if (gt1Post < gt0Post) {
				  dest.setGenotypeQuality(sample, gt1Post - gt2Post);
				  return 2;
				}
				else {
				  dest.setGenotypeQuality(sample, gt0Post - gt2Post);
				  return 2;
				}
			}
			else {
				dest.setGenotypeQuality(sample, gt2Post - gt0Post);
			    return 0;
			}
		}
		else if (gt1Post < gt0Post) {
			if (gt0Post < gt2Post) {
				dest.setGenotypeQuality(sample, gt0Post - gt1Post);
			    return 1;
			}
			else {
			  dest.setGenotypeQuality(sample, gt2Post - gt1Post);
			  return 1;
			}
		}
		else {
			dest.setGenotypeQuality(sample, gt1Post - gt0Post);
		    return 0;
		}
	}

	private final static double C1P = - 10.0 / Math.log(10);

	private static double phredOf(double p) {
		if (p > 0.99) 
			return C1P * Math.log1p(p - 1);
		else
			return - 10 * Math.log10(p);
	}
	
	private static double phredOf1p(double p) {
		if (p < 0.01)
			return C1P * Math.log1p(- p);
		else
			return - 10 * Math.log10(1 - p);
	}
	
	
	private class PriorHeader {
		private PriorFormat format;
		private int alleleCount;
		private int categoryCount;
		private PriorEncoding encoding;
		private double base = -1;
		private double scale = -1;
		
	
		public PriorHeader(String line) {
			if (line == null)
				throw new GenotypingModelException(String.format("no header found in prior file '%s'",prior.toString()));
			String[] headerParts = line.split("\\s+");
			if (headerParts.length == 0)
				throw new GenotypingModelException(String.format("empty header found in prior file '%s'",prior.toString()));
			if (headerParts.length < 3)
				throw new GenotypingModelException(String.format("incomplete header found in prior file '%s'", prior.toString()));
			if (! headerParts[0].startsWith("#"))
				throw new GenotypingModelException(String.format("header does no start with '#' in prior file '%s'", prior.toString()));
			if ((format = PriorFormat.valueOf(headerParts[0].substring(1))) == null)
				throw new GenotypingModelException(String.format("unknown prior file format '%s' found in file '%s'",headerParts[0],prior.toString()));
			alleleCount = readHeaderInteger(headerParts[1],"allele count");
			if (alleleCount < 2)
				throw new GenotypingModelException(String.format("prior allele count must be at least 2"));
			categoryCount = readHeaderInteger(headerParts[2],"category count");
			if (categoryCount < 2)
				throw new GenotypingModelException(String.format("prior category count must be at least 2"));
			if ((encoding = headerParts.length < 4 ? PriorEncoding.LINEAR : PriorEncoding.valueOf(headerParts[3]) ) == null)
				throw new GenotypingModelException(String.format("illegal encoding '%s'",headerParts[3]));
			if (encoding == PriorEncoding.PHRED && headerParts.length > 4)
                throw new GenotypingModelException(String.format("Phred encoding does not accept further paraters (base and scale)"));
			if (headerParts.length > 4) {
				base = readHeaderDouble(headerParts[4],"base");
				if (base <= 0 || Double.isInfinite(base) || Double.isNaN(base))
					throw new GenotypingModelException(String.format("illegal base '%s'",headerParts[4]));
			}
			else {
				switch (encoding) {
				case LINEAR:
					base = 0; break;
 				case PHRED:
 					base = 10; break;
 				case LN:
 					base = Math.E; break;
 				case LOG:
 					throw new GenotypingModelException("with encoding LOG you need to specify a base");
				}
			}
			if (headerParts.length > 5) {
				scale = readHeaderDouble(headerParts[5],"scale");
				if (Double.isInfinite(scale) || Double.isNaN(scale))
					throw new GenotypingModelException(String.format("illegal scale '%s'",headerParts[5]));
			}
			else {
				scale = encoding == PriorEncoding.PHRED ? -10 : 1;
			}
				
			
		}
		
		public double readLine(String line, double[] priorProbs,
				double[] priorFreqs, double[] oneMinusPriorProbs, int idx) {
			String[] parts = line.split("\\s+");
			if (parts.length != 2)
				throw new GenotypingModelException("illegal prior file content at category '" + idx + "'");
			double freq = readDouble(parts[0],String.format("category %d frequency",idx+1));
			if (Double.isNaN(freq))
				throw new GenotypingModelException("NaN not a allowed as a frequency in prior");
			if (freq < 0)
				throw new GenotypingModelException("frequencies cannot be negative in prior");
			if (freq > 1)
				throw new GenotypingModelException("frequencies cannot be greater than one in prior");
		    priorFreqs[idx] = freq;
		    
		    double probs = readDouble(parts[1],String.format("category %d probability",idx+1));
		    if (Double.isNaN(probs))
		    	throw new GenotypingModelException("NaN not allowed as a probability in priors");
		    if (Double.isInfinite(probs))
		    	throw new GenotypingModelException("Infinite not allowedd as a probability in priors");    
		    if (probs < 0)
		    	throw new GenotypingModelException("Probability");
		    double oneMinusProbs = Double.NaN;
		    double probsDivScale = scale == 1.0 ? probs : probs/scale;
		    switch (this.encoding) {
		    case LINEAR:
		    	probs = phredOf(probsDivScale);
		    	oneMinusProbs = phredOf1p(probsDivScale);
		    	break;
		    case PHRED:
		    	oneMinusProbs = phredOf1p(Math.pow(10,-probs/10));
		    	break;
		    case LN:
		    	probs = phredOf(Math.exp(probsDivScale));
		    	oneMinusProbs = phredOf1p(Math.exp(probsDivScale));
		    	break;
		    case LOG:
		    	probs = phredOf(Math.pow(base,probsDivScale));
		    	oneMinusProbs = phredOf1p(Math.pow(base, probsDivScale));
		    }
		    priorProbs[idx] = probs;
		    oneMinusPriorProbs[idx] = oneMinusProbs;
		    
		    return Math.pow(10,-probs/10);
		}
		

		private int readHeaderInteger(String str,String role) {
			try {
				return Integer.parseInt(str);
			}
			catch (NumberFormatException e) {
				throw new GenotypingModelException(String.format("prior header in file '%s' contains an invalid %s (%s)",prior.toString(),role,str));
			}
		}

		private double readDouble(String str,String role) {
			try {
				if (str.equalsIgnoreCase("e")) {
					return Math.E;
				}
				else if (str.equalsIgnoreCase("pi")) {
					return Math.PI;
				}
				return Double.parseDouble(str);
			}
			catch (NumberFormatException e) {
				throw new GenotypingModelException(String.format("prior in file '%s' contains an invalid %s (%s)",prior.toString(),role,str));
			}
		}		
		
		private double readHeaderDouble(String str,String role) {
			try {
				if (str.equalsIgnoreCase("e")) {
					return Math.E;
				}
				else if (str.equalsIgnoreCase("pi")) {
					return Math.PI;
				}
				return Double.parseDouble(str);
			}
			catch (NumberFormatException e) {
				throw new GenotypingModelException(String.format("prior header in file '%s' does not contain a valid %s (%s)",prior.toString(),role,str));
			}
		}

	}
	
	private enum PriorFormat {
		PRIOR	
	}
	
	private enum PriorEncoding {
		LINEAR, LOG, LN, PHRED
	}
	
	
	public static final String PFREQ_AVERAGE_KEY = "PF";
	public static final String PFREQ_SD_KEY = "SF";	
    public static final String BAYES_FACTORS_KEY = "BF";
    
    private static final VCFFormatHeaderLine PFREQ_FORMAT_LINE = 
    		new VCFFormatHeaderLine(PFREQ_AVERAGE_KEY, 1, VCFHeaderLineType.Float, "Average posterior alternative allele frequency [0..1]");
    
    private static final VCFFormatHeaderLine PFREQ_SD_FORMAT_LINE = 
    		new VCFFormatHeaderLine(PFREQ_SD_KEY, 1, VCFHeaderLineType.Float, "Posterior alternaive allele frequency std. deviation.");
	
	private static final VCFFormatHeaderLine BAYES_FACTORS_FORMAT_LINE = 
			new VCFFormatHeaderLine(BAYES_FACTORS_KEY, -1, VCFHeaderLineType.Float, "Bayes factors for possible genotypes #0 vs #1, #0 vs #2 ... #0 vs #N, #1 vs #2 ... #N-1 vs #N");
    
    @Override
    public Collection<? extends VCFHeaderLine> getHeaderLines() {
		List<VCFHeaderLine> list = new LinkedList<VCFHeaderLine>();
		list.addAll(super.getHeaderLines());
		list.add(PFREQ_FORMAT_LINE);
		list.add(PFREQ_SD_FORMAT_LINE);
		list.add(BAYES_FACTORS_FORMAT_LINE);
		return Collections.unmodifiableList(list);
    }
}
