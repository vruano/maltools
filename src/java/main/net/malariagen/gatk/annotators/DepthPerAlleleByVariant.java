package net.malariagen.gatk.annotators;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class DepthPerAlleleByVariant extends InfoFieldAnnotation {

	public static final List<String> KEY_NAMES = Collections
			.singletonList("AD");
	public static final List<VCFInfoHeaderLine> DESCRIPTIONS = Collections
			.singletonList(new VCFInfoHeaderLine(KEY_NAMES.get(0),
					-1, VCFHeaderLineType.Integer,
					"Depth per allele across all samples"));

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			AnnotatorCompatibleWalker walker, ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {

		if (!vc.isSNP())
			return null;
		byte[] bases = new byte[vc.getAlleles().size()];
		if (bases.length == 0)
			return null;
		bases[0] = vc.getReference().getBases()[0];
		int[] counts = new int[bases.length];
		for (int i = 1; i < bases.length; i++)
			bases[i] = vc.getAlternateAllele(i - 1).getBases()[0];
		for (AlignmentContext ac : stratifiedContexts.values())
			for (PileupElement pe : ac.getBasePileup()) {
				byte b = pe.getBase();
				for (int i = 0; i < bases.length; i++)
					if (bases[i] == b) {
						counts[i]++;
						break;
					}
			}
		return Collections.singletonMap(KEY_NAMES.get(0), (Object) counts);
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
