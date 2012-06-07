package net.malariagen.gatk.annotators;

import java.util.Collections;
import java.util.List;
import java.util.Map;


import net.malariagen.gatk.uniqueness.UQNFeature;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class UniquenessScore implements InfoFieldAnnotation {

	private static final Integer NO_SCORE = 99;
	
	// cached Integer object representing the uniqueness scores.
	private static final Integer[] SCORES;
	
	static {
		SCORES = new Integer[NO_SCORE + 1];
		for (int i = 0; i <= NO_SCORE; i++)
			SCORES[i] = i;
	}

	private static List<String> KEY_NAMES = Collections.singletonList(Constants.UNIQUENESS_KEY);
	
	private static List<VCFInfoHeaderLine> HEADER_LINES = Collections.singletonList(
			new VCFInfoHeaderLine(Constants.UNIQUENESS_KEY, 1, VCFHeaderLineType.Integer, 
					"Indicates the uniqueness of the genomic region where the variant is found"));
	
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		Integer score = score(tracker);
		if (score == null)
			return Collections.emptyMap();
		else
		  return Collections.singletonMap(Constants.UNIQUENESS_KEY, (Object) score );		
	}
	
	private Integer score(RefMetaDataTracker tracker) {
		List<GATKFeature> tracks = tracker.getGATKFeatureMetaData(Constants.UNIQUENESS_ROD_NAME, true);
		for (GATKFeature ft : tracks) {
			Object o = ft.getUnderlyingObject();
			if (o instanceof UQNFeature) {
				UQNFeature u = (UQNFeature)o;
				return SCORES[u.getScore()];
			}
		}
		return NO_SCORE;
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
