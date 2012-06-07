package net.malariagen.gatk.annotators;

import static net.malariagen.gatk.annotators.Constants.*;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class AbsoluteCounts extends AlleleCountsAnnotator implements InfoFieldAnnotation {
	

	private List<String> KEY_NAMES = Collections.unmodifiableList(Arrays.asList(ABSOLUTE_ALLELE_DEPTH_KEY,ABSOLUTE_TOTAL_DEPTH_KEY,MINOR_READ_ALLELE_KEY,MINOR_READ_ALLELE_FREQUENCY_KEY,"mAmaxSF"));

	private List<VCFInfoHeaderLine> HEADER_LINES = Collections.unmodifiableList(Arrays.asList(
	    new VCFInfoHeaderLine(ABSOLUTE_ALLELE_DEPTH_KEY, VCFInfoHeaderLine.UNBOUNDED, VCFHeaderLineType.Integer, "Non-filtered indel-free depth per allele"),
	    new VCFInfoHeaderLine(ABSOLUTE_ALLELE_MAX_SAMPLE_DEPTH_KEY, VCFInfoHeaderLine.UNBOUNDED, VCFHeaderLineType.Integer, "Non-filtered indel-free depth maximum per sample depth"),
		new VCFInfoHeaderLine(ABSOLUTE_TOTAL_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Non-filtered indel-free depth per allele"),
		new VCFInfoHeaderLine(MINOR_READ_ALLELE_KEY, 1, VCFHeaderLineType.String, "Minor read Allele (considering absolute read depth)"),
		new VCFInfoHeaderLine(MINOR_READ_ALLELE_TOTAL_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Minor read Allele depth of coverage"),
		new VCFInfoHeaderLine(MINOR_READ_ALLELE_FREQUENCY_KEY, 1, VCFHeaderLineType.String, "Minor read Allele Frequency (as a fraction [0 to 0.5))"),
		new VCFInfoHeaderLine(LOWEST_MAXIMUM_READ_SAMPLE_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Lowest maximum allele sample Depth: the lowest maximum depth of any particular allele in any sample"),
		new VCFInfoHeaderLine(MINOR_READ_ALLELE_MAXIMUM_SAMPLE_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Minor read Allele maximum per Sample Depth: absolute depth of the MrA at the sample in which it is more frequently observed")));
	
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		
		int[] alleleCounts;
		int[] maxSampleDepths;
		
		char[] ac = initializeAlleleCharsArray();
		int[] ai = initializeAlleleIndecesArray();
		int alleleCount = initializeAlleleArrays(vc, ac, ai);

		alleleCounts = new int[alleleCount];
		maxSampleDepths = new int[alleleCount];
		int[] sampleDepths = new int[alleleCount];
		
        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
        	AlignmentContext alignmentContext = sample.getValue();
        	Arrays.fill(sampleDepths, 0);
        	if ( ! alignmentContext.hasBasePileup() ) continue;
	        ReadBackedPileup pileup = alignmentContext.getBasePileup();
	        for ( PileupElement p : pileup ) {
	        	int b = p.getBase();
	        	int i = ai[b - Byte.MIN_VALUE];
	        	if (i == -1) continue;
	        	sampleDepths[i] += 1;
	        }
	        for (int i = 0; i < alleleCount; i++) {
	        	int sd = sampleDepths[i];
	        	alleleCounts[i] += sd;
	        	if (sd > maxSampleDepths[i]) maxSampleDepths[i] = sd;
	        }
        }
        int minorReadAlleleIndex = 0;
        int lowestMaximumReadSampleDepthIndex = 0;
        int total = alleleCounts[0];
        for (int i = 1; i < alleleCount; i++) {
        	total += alleleCounts[i];
        	if (alleleCounts[i] < alleleCounts[minorReadAlleleIndex])
        		minorReadAlleleIndex = i;
        	if (maxSampleDepths[i] < maxSampleDepths[lowestMaximumReadSampleDepthIndex]) 
        		lowestMaximumReadSampleDepthIndex = i;
        }

       
        Map<String, Object> map = new HashMap<String,Object>(2);
        map.put(ABSOLUTE_ALLELE_DEPTH_KEY, alleleCounts);
        map.put(ABSOLUTE_TOTAL_DEPTH_KEY, total);
        map.put(LOWEST_MAXIMUM_READ_SAMPLE_DEPTH_KEY, alleleCounts[lowestMaximumReadSampleDepthIndex]);
        map.put(ABSOLUTE_ALLELE_MAX_SAMPLE_DEPTH_KEY, maxSampleDepths);
        map.put(MINOR_READ_ALLELE_KEY, "" + ac[minorReadAlleleIndex]);
        map.put(MINOR_READ_ALLELE_TOTAL_DEPTH_KEY, "" + alleleCounts[minorReadAlleleIndex]);
        map.put(MINOR_READ_ALLELE_MAXIMUM_SAMPLE_DEPTH_KEY, maxSampleDepths[minorReadAlleleIndex]);
        double mAF = ((double)alleleCounts[minorReadAlleleIndex]) / ((double)total);
        map.put(MINOR_READ_ALLELE_FREQUENCY_KEY, formatFractionAnnotation(mAF));
        return map;
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
