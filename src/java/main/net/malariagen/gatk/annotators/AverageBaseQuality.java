package net.malariagen.gatk.annotators;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.util.FastMath;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class AverageBaseQuality implements InfoFieldAnnotation {
	public static final List<String> KEY_NAMES = Collections
			.singletonList("ER");
	public static final List<VCFInfoHeaderLine> DESCRIPTIONS = Collections
			.singletonList(new VCFInfoHeaderLine(KEY_NAMES.get(0),
					1, VCFHeaderLineType.Float,
					"E[min(BQ,MQ)] i.e. the average base quality considering the minimum between base call and including read mapping quality"));

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		int size = 0;
		for (AlignmentContext ac : stratifiedContexts.values())
			size += ac.size();
				
		if (size == 0)
			return null;
				
		int[] quals = new int[size];
		int nextIdx = 0;
		
		for (AlignmentContext ac : stratifiedContexts.values()) 
			for (PileupElement pe : ac.getBasePileup()) {
				byte bq = pe.getQual();
				int mq = pe.getMappingQual();
				int aq = bq > mq ? mq : bq;
				quals[nextIdx++] = aq;
				FastMath.min(pe.getQual(),pe.getMappingQual());
			}
		
		double expSum = 0;
		for (int i = 0; i < quals.length; i++) {
			expSum += FastMath.pow(10, ((double) - quals[i]) / 10.0);
		}
		
				
		//double result = -10 * FastMath.log10(expSum) - c + 10 * FastMath.log10(quals.length);
		double result = -10 * FastMath.log10(expSum/(double)quals.length);
			
		String resultString = String.format("%.2f", result);
				
		return Collections.singletonMap(KEY_NAMES.get(0), (Object) resultString);
	}

	@Override
	public List<String> getKeyNames() {
		return KEY_NAMES;
	}

	@Override
	public List<VCFInfoHeaderLine> getDescriptions() {
		return DESCRIPTIONS;
	}
}
