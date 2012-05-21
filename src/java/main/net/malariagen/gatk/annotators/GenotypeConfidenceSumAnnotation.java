package net.malariagen.gatk.annotators;


import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class GenotypeConfidenceSumAnnotation implements InfoFieldAnnotation {
	

	public static final String GENOTYPE_CONFIDENCE_SUM_KEY = "GCSum";
	
	private static List<String> KEY_NAMES = Collections.unmodifiableList(Arrays.asList(GENOTYPE_CONFIDENCE_SUM_KEY));

	private static double MIN_CONFIDENCE = -99999;
	private static List<VCFInfoHeaderLine> HEADER_LINES = Collections.singletonList(
		new VCFInfoHeaderLine(GENOTYPE_CONFIDENCE_SUM_KEY, -1, VCFHeaderLineType.Float, "Genotype confidence sum"));
	
	
	@Override
	public Map<String,Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> acs, VariantContext vc) {

		double sum = 0;
		for (Genotype g : vc.getGenotypes().values()) {
			Object o = g.getAttribute("GC");
			if (o instanceof CharSequence)
				sum += Double.parseDouble(o.toString());
			else if (o instanceof Number)
				sum += ((Number)o).doubleValue();
		}
		if (sum < MIN_CONFIDENCE)
			sum = MIN_CONFIDENCE;
		return Collections.singletonMap(GENOTYPE_CONFIDENCE_SUM_KEY, (Object) new Double(sum));
	}


	@Override
	public List<String> getKeyNames() {
		return KEY_NAMES;
	}

	@Override
	public List<VCFInfoHeaderLine> getDescriptions() {
		return HEADER_LINES;
	}

}
