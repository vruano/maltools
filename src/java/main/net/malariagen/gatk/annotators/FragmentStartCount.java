package net.malariagen.gatk.annotators;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class FragmentStartCount extends InfoFieldAnnotation {
	
	private static final String FORWARD_START_KEY = "FFS";
	private static final String REVERSE_START_KEY = "FRS";
	private static final List<String> KEY_NAMES = Arrays.asList(FORWARD_START_KEY,REVERSE_START_KEY);
	private static final List<VCFInfoHeaderLine> DESCRIPTIONS = Arrays.asList(
			new VCFInfoHeaderLine(FORWARD_START_KEY,1,VCFHeaderLineType.Integer,"Number of fragment forward starts found"),
			new VCFInfoHeaderLine(REVERSE_START_KEY,1,VCFHeaderLineType.Integer,"Number of fragmen reverse starts found")
	);
	private final int minMappingQuality;
	
	public FragmentStartCount() {
		this(0);
	}
	
	public FragmentStartCount(int minMQ) {
		minMappingQuality = minMQ;
	}
	
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			AnnotatorCompatibleWalker walker, ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		
		int forwardStarts = 0;
		int reverseStarts = 0;
		for (AlignmentContext ac : stratifiedContexts.values())
			for (PileupElement pe : ac.getBasePileup()) {
				if (pe.getMappingQual() < minMappingQuality)
					continue;
				SAMRecord read = pe.getRead();
				if (!read.getFirstOfPairFlag()) continue;
				if (pe.getOffset() == 0) {
					if (!read.getReadNegativeStrandFlag()) forwardStarts++;
				}
				else if (pe.getOffset() == pe.getRead().getReadLength() - 1) {
					if (read.getReadNegativeStrandFlag()) reverseStarts++;
				}
			
			}
		Map<String,Object> result = new HashMap<String,Object>(2);
		result.put(FORWARD_START_KEY,forwardStarts);
		result.put(REVERSE_START_KEY,reverseStarts);
		return result;
	}

	@Override
	public List<VCFInfoHeaderLine> getDescriptions() {
		return DESCRIPTIONS;
	}

	@Override
	public List<String> getKeyNames() {
		return KEY_NAMES;
	}

}
