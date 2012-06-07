package net.malariagen.gatk.genotyper;

import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;

public class GenotypePosteriors {

	protected double[] posteriors;

	public GenotypePosteriors(double[] p) {
		posteriors = p;
	}

	@Override
	public String toString() {
		if (posteriors == null || posteriors.length == 0)
			return VCFConstants.MISSING_VALUE_v4;
		StringBuffer sb = new StringBuffer();
		for (double d : posteriors) {
			if (d == 0.0) d = 0.0;
			sb.append(String.format(VCFConstants.DOUBLE_PRECISION_FORMAT_STRING, d)).append(',');
		}
		sb.setLength(sb.length() - 1);
		return sb.toString();
	}

	public double getPosterior(int i) {
		return posteriors[i];
	}

}
