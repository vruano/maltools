package net.malariagen.gatk.annotators;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class NumberSamplesWithData implements InfoFieldAnnotation {
	public static final List<String> KEY_NAMES = Collections
			.singletonList("NS");
	public static final List<VCFInfoHeaderLine> DESCRIPTIONS = Collections
			.singletonList(new VCFInfoHeaderLine(KEY_NAMES.get(0),
					1, VCFHeaderLineType.Integer,
					"Number of samples with data after all filtering"));

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		int coveredSamples = 0;
		for (AlignmentContext ac : stratifiedContexts.values())
			if (ac.size() > 0) coveredSamples++;
		return Collections.singletonMap(KEY_NAMES.get(0), (Object) coveredSamples);
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
