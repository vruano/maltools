package net.malariagen.gatk.genotyper.models;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.genotyper.GenotypePosteriors;
import net.malariagen.gatk.genotyper.VariantPosteriors;
import net.malariagen.gatk.genotyper.GenotypingContext;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.InferredGeneticContext;
import org.broadinstitute.sting.utils.variantcontext.MutableGenotype;


public abstract class AbstractGenotypingModel implements GenotypingModel {
	private GenotypingContext gc;
	
	public static final String ERROR_RATE_KEY = "ER";
	public static final String GENOTYPE_CONFIDENCE_KEY = "GC";

	public static final VCFFormatHeaderLine ERROR_RATE_FORMAT_LINE = new VCFFormatHeaderLine(ERROR_RATE_KEY, 1, VCFHeaderLineType.Float, "Average base error/quality of considered calls (cap by read mapping quality and BAQ if applied)");
	public static final VCFFormatHeaderLine GENOTYPE_CONFIDENCE_FORMAT_LINE = new VCFFormatHeaderLine(GENOTYPE_CONFIDENCE_KEY, 1, VCFHeaderLineType.Float, "Genotype call confidence: indicates how well the data evidence fits the model assumed");
    public static final VCFFormatHeaderLine GENOTYPE_POSTERIORS_FORMAT_LINE = new VCFFormatHeaderLine(VCFConstants.GENOTYPE_POSTERIORS_KEY, -1, VCFHeaderLineType.Float, "Each possible genotype posterior probability, #0, #1, #2 ... #N");
	
	
	private static final Collection<? extends VCFHeaderLine> HEADER_LINES = Collections.unmodifiableCollection(
			Arrays.asList(ERROR_RATE_FORMAT_LINE,GENOTYPE_CONFIDENCE_FORMAT_LINE, GENOTYPE_POSTERIORS_FORMAT_LINE)); 
	
	@Override
	public Collection<? extends VCFHeaderLine> getHeaderLines() {
		return HEADER_LINES;
	}
	
	public void setGenotypingContext(GenotypingContext gc) {
		if (gc.equals(this.gc))
			return;
		this.gc = gc;
		init();
	}
	
	@Override
	public String getModelName() {
		GModel anno = this.getClass().getAnnotation(GModel.class);
		if (anno != null) { 
			if (!anno.shortName().equals(GModel.NO_NAME))
				return anno.shortName();
			else if (!anno.fullName().equals(GModel.NO_NAME))
				return anno.fullName();
			else if (anno.aliases().length > 0) 
				return anno.aliases()[0];
		}
		return this.getClass().getSimpleName();
	}
	
	void init() {
		GenotypingContext gc = this.getGenotypingContext();
		if (gc == null)
			throw new IllegalStateException("no genotyping context has been provided to te model");
	}
	
	public GenotypingContext getGenotypingContext() {
		return gc;
	}
	
	protected VariantPosteriors calculatePosteriors(AlignmentContext ... acs) {
		VariantPosteriors result = new VariantPosteriors(acs.length,getGenotypeCount());
		for (int i = 0; i < acs.length; i++) {
			AlignmentContext ac = acs[i];
			calculatePosteriors(ac, result, i);
		}
		return result;
	}

	protected abstract double calculatePosteriors(AlignmentContext ac, VariantPosteriors dest, int sample);

	@Override
	public Map<String, MutableGenotype> calculateGenotypes(Map<String, AlignmentContext> sc) {
		String[] samples = sc.keySet().toArray(new String [sc.size()]);
		AlignmentContext[] acs = new AlignmentContext[samples.length];
		for (int i = 0; i < samples.length; i++)
			acs[i] = sc.get(samples[i]);
		VariantPosteriors posteriors = calculatePosteriors(acs);
		HashMap<String, MutableGenotype> result = new HashMap<String,MutableGenotype>(acs.length);
		for (int i = 0; i < acs.length; i++)  
			result.put(samples[i],calculateGenotype(samples[i],acs[i],posteriors,i));
		return result;
	}

	public MutableGenotype calculateGenotype(String sample,
			AlignmentContext ac, VariantPosteriors post, int i) {
		double genotypeQuality =  post.getGenotypeQuality(i);
		if (genotypeQuality > 99999) { genotypeQuality = 99999;  }
		else if (genotypeQuality == -0.0) { genotypeQuality = 0.0; }
		double negLog10PError = Double.isNaN(genotypeQuality) ? InferredGeneticContext.NO_NEG_LOG_10PERROR : genotypeQuality / 10.0;
		return new MutableGenotype(sample, genotypeAlleles(ac, post, i),
				negLog10PError,
				genotypeFilters(ac, post, i), genotypeAttributes(ac,
						post, i), false);
	}

	public Map<String, ?> genotypeAttributes(AlignmentContext ac, VariantPosteriors post, int i) {
		Map<String, Object> result = new HashMap<String, Object>();
		StringBuffer sb = new StringBuffer();
		post.getPosteriorsString(i, sb);
		result.put(VCFConstants.GENOTYPE_POSTERIORS_KEY, post.getGenotypePosteriors(i));
		if (!Double.isNaN(post.getGenotypeConfidence(i)))
				result.put("GC",String.format("%.2f",post.getGenotypeConfidence(i)));
		if (!Double.isNaN(post.getAvgErrRate(i)))
			result.put("ER",String.format("%.2f",post.getAvgErrRate(i)));
		//Double pf = (Double) post.getAttribute(i,"PF");
		//if (pf != null)
		//	result.put("PF", pf > 0.00005 ? "0" : String.format("%.2g",pf));
		//Double sf = (Double) post.getAttribute(i,"SF");
		//if (sf != null)  result.put("SF", sf < 0.00005 ? "0" : String.format("%.2g",sf));
		
		String bf = (String) post.getAttribute(i,"BF");
		if (bf != null) result.put("BF", bf);
		return result;
	}

	public Set<String> genotypeFilters(AlignmentContext ac, VariantPosteriors post, int i) {
		return null;
	}

	public List<Allele> genotypeAlleles(AlignmentContext ac, VariantPosteriors post, int i) {
		int genotypeCount = getGenotypeCount();
		int bestGenotype = -1;
		double bestScore = Integer.MAX_VALUE;
		for (int j = 0; j < genotypeCount; j++) {
			double score = post.getPosterior(i, j);
			if (score < 0) continue;
			if (score < bestScore) {
				bestGenotype = j;
				bestScore = score;
			}
		}
		return getGenotypeAlleles(bestGenotype);
	}
	
	@Override
	public double calculateVariantPhredQuality(Map<String, AlignmentContext> ac, Map<String, Genotype> gt) {
		double value = 0;
		for (Genotype g : gt.values()) {
			if (!g.isCalled()) continue;
			int idx = this.getGenotypeIndex(g.getAlleles());
			if (idx == -1)
				continue;
			GenotypePosteriors gp = (GenotypePosteriors) g.getAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY);
			if (gp == null)
				continue;
			double v = gp.getPosterior(0);
			if (!Double.isNaN(v) && v != InferredGeneticContext.NO_NEG_LOG_10PERROR)
				value += v;
		}
		if (value > 99999.99) value = 99999.99;
		return value;
	}

}
