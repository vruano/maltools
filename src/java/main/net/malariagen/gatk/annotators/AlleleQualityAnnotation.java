package net.malariagen.gatk.annotators;


import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.utils.NucleotideIUPAC;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class AlleleQualityAnnotation implements InfoFieldAnnotation {
	
	public static final String ALLELE_QUALITY_KEY = "AQ";
	
	private static List<String> KEY_NAMES = Collections.unmodifiableList(Arrays.asList(ALLELE_QUALITY_KEY));


	private static double MAX_QUALITY = 9999;
	private static List<VCFInfoHeaderLine> HEADER_LINES = Collections.singletonList(
		new VCFInfoHeaderLine(ALLELE_QUALITY_KEY, -1, VCFHeaderLineType.Float, "Believe in that the allele actually exist, as long as mapping is correct"));
	
	@Override
	public Map<String,Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> acs, VariantContext vc) {
		
		if (!vc.isSNP())
			return Collections.emptyMap();
		Set<Allele> alleles = vc.getAlternateAlleles();
		
		if (alleles.size() == 0)
			return Collections.emptyMap();
		NucleotideCounts ncs = new NucleotideCounts();
		ncs.add(acs.values());
		ncs.updateQualities();
		StringBuffer sb = new StringBuffer();
		for (Allele a : alleles) {
			double quality = ncs.getQuality(a);
			if (quality > MAX_QUALITY) quality = MAX_QUALITY;
			sb.append(String.format("%.2f",quality)).append(',');	
		}
		if (sb.length() > 0) 
			sb.setLength(sb.length() -1);
		
		return Collections.singletonMap(ALLELE_QUALITY_KEY, (Object) sb.toString());
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
