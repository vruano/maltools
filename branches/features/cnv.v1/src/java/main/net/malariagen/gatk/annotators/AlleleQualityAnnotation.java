package net.malariagen.gatk.annotators;


import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.math.Beta;
import net.malariagen.utils.NucleotideIUPAC;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
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
		double errorRateSum = 0;
		int[] alleleCounts = new int[NucleotideIUPAC.values().length];
		int total = 0;
		for (AlignmentContext ac : acs.values()) {
			ReadBackedPileup rb = ac.getBasePileup();
			for (PileupElement e : rb) {
				alleleCounts[NucleotideIUPAC.fromBase(e.getBase()).ordinal()]++;
				byte qual = baseQual(e);
				total++;
				errorRateSum += Math.pow(10, -qual * 0.1);
			}
		}
		double errorRateAvg = errorRateSum / total;
		StringBuffer sb = new StringBuffer(10 * alleles.size());
		for (Allele a : alleles) {
			int count = alleleCounts[NucleotideIUPAC.fromBase(a.getBases()[0]).ordinal()];
			double quality = alleleQuality(errorRateAvg,total,count);
			if (quality > MAX_QUALITY)
				quality = MAX_QUALITY;
			sb.append(String.format("%.2f",quality)).append(',');	
		}
		if (sb.length() > 0) 
			sb.setLength(sb.length() -1);
		return Collections.singletonMap(ALLELE_QUALITY_KEY, (Object) sb.toString());

	}
	
	public double alleleQuality (double er, int total, int count) {	
	    if (count == 0) 
	       return - 10 * Beta.log10(er,1.0000001,total + 0.000001);
	    else if (total == count)
	       return - 10 * Beta.log10(1 - er, 1.000001,count + 0.000001);
	    else
	       return - 10 * Beta.log10(er,count + 1.000001,total - count + 0.000001);
	    
	}
	
	private byte baseQual(PileupElement e) {
		int mq = e.getMappingQual();
		byte bq = e.getQual();
		byte result = mq < bq ? (byte) mq : bq;
		// 3 equal to 0.5 chances. 
		return result < 3 ? 3 : result;
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
