package net.malariagen.gatk.annotators;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import net.malariagen.gatk.gff.GFFFeature;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class CodingAnnotation implements InfoFieldAnnotation {

	private static List<String> KEY_NAMES = Collections.singletonList(Constants.CODING_KEY);
	
	private static List<VCFInfoHeaderLine> HEADER_LINES = Collections.singletonList(
			new VCFInfoHeaderLine(Constants.CODING_KEY, 0, VCFHeaderLineType.Flag, "Mark variants at protein coding positions"));
	
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		return Collections.singletonMap(Constants.CODING_KEY,(Object)isCoding(tracker));		
	}
	
	static boolean isCoding(RefMetaDataTracker tracker) {
		List<GATKFeature> tracks = tracker.getGATKFeatureMetaData(Constants.FEATURES_ROD_NAME, true);
		for (GATKFeature ft : tracks) {
			Object o = ft.getUnderlyingObject();
			if (o instanceof GFFFeature) {
				if (((GFFFeature)o).getType().isProteinCoding()) return true;
			}
		}
		return false;
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
