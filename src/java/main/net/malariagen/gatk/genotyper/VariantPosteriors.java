package net.malariagen.gatk.genotyper;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.util.FastMath;

public class VariantPosteriors {
	
	protected double[] values;
	int genotypeCount;
	int sampleCount;
	
	Map<Integer, Map<String, Object>> attributes; 
	
	
	protected double[] qualities;
	protected double[] confidence;
	protected double[] avgErrRate;

	public VariantPosteriors (int sampleCount, int genotypeCount) {
		checkDimensions(sampleCount,genotypeCount);
		values = new double[sampleCount * genotypeCount];
		Arrays.fill(values,Double.NaN);
		resize(sampleCount,genotypeCount);
		qualities = new double[sampleCount];
		confidence = new double[sampleCount];
		avgErrRate = new double[sampleCount];
		Arrays.fill(qualities, Double.NaN);
		Arrays.fill(confidence, Double.NaN);
		Arrays.fill(avgErrRate, Double.NaN);
		clear();
	}
	
	public void setAttribute(int sample, String name, Object value) {
		if (attributes == null)
			attributes = new HashMap<Integer,Map<String,Object>>(sampleCount);
		Map<String,Object> sampleAttributes = attributes.get(sample);
		if (sampleAttributes == null)
			attributes.put(sample, sampleAttributes = new HashMap<String,Object>(10));
	    sampleAttributes.put(name,value);
	}
	
	public Object getAttribute(int sample, String name) {
		if (attributes == null)
			return null;
		Map<String,Object> sampleAttributes = attributes.get(sample);
		if (sampleAttributes == null)
			return null;
		return sampleAttributes.get(name);
	}
	
	
	public void setAvgErrRate(int sample, double aer) {
		avgErrRate[sample] = aer;
	}
	
	public double getAvgErrRate(int sample) {
		return avgErrRate[sample];
	}
	
	
	private void checkDimensions(int sampleCount, int genotypeCount) {
		if (sampleCount < 0 || genotypeCount < 0)
			throw new IllegalArgumentException("Illegal Dimensions");
	}

	public void clear() {
		Arrays.fill(values, Integer.MAX_VALUE);
	}
	
	public void resize(int sampleCount, int genotypeCount) {
		this.sampleCount = sampleCount;
		this.genotypeCount = genotypeCount;
		if (sampleCount * genotypeCount > values.length) {
			values = new double[sampleCount * genotypeCount];
			clear();
		}
	}
	
	public double getGenotypeQuality(int sample) {
		return qualities[sample];
	}
	
	public void setGenotypeConfidence(int sample, double conf) {
		confidence[sample] = conf;
	}
	
	public double getGenotypeConfidence(int sample) {
		return confidence[sample];
	}
	
	public void setGenotypeQuality(int sample, double q) {
		qualities[sample] = q;
	}
	
	public void setGenotypePosterior(int sample, int genotype, double value) {
		values[sample * genotypeCount + genotype] = value;
	}
	
	public double getPosterior(int sample, int genotype) {
		return values[sample * genotypeCount + genotype];
	}
	
	@Deprecated
	public void getPosteriorsString(int sample, StringBuffer sb) {
		sb.append(getGenotypePosteriors(sample).toString());
	}
	
	public GenotypePosteriors getGenotypePosteriors(int sample) {
		double[] p = new double[genotypeCount];
		System.arraycopy(values, sample * genotypeCount, p, 0, genotypeCount);
		return new GenotypePosteriors(p);
	}
	
	public String getSamplePosteriorsString(int sample) {
		StringBuffer sb = new StringBuffer();
		getPosteriorsString(sample,sb);
		return sb.toString();
	}
	
	public void getSamplePosteriors(int sample, double[] dest, int offset) {
		if (dest.length - offset < genotypeCount)
			throw new IllegalArgumentException("there is no enough space in destination array at indicated offset");
		System.arraycopy(values,sample * genotypeCount, dest, offset, genotypeCount);
	}
	
	public double[] getSamplePosteriors(int sample) {
		double[] result = new double[genotypeCount];
		getSamplePosteriors(sample,result,0);
		return result;
	}
	
	public void getSamplePosteriors(int sample, double[] dest) {
		getSamplePosteriors(sample,dest,0);
	}
	
	public void setSamplePosteriors(int sample, double ... values) {
		System.arraycopy(values, 0, this.values, sample * genotypeCount, values.length);
	}

}
