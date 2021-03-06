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
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class MinorAlleleCounts extends  InfoFieldAnnotation {
	
	
	private static List<String> KEY_NAMES = Collections.unmodifiableList(Arrays.asList(MINOR_ALLELE_KEY,MINOR_ALLELE_FREQUENCY_KEY));

	private static List<VCFInfoHeaderLine> HEADER_LINES = Collections.unmodifiableList(Arrays.asList(
		new VCFInfoHeaderLine(MINOR_ALLELE_KEY, 1, VCFHeaderLineType.String, "Minor Allele based in individual sample genotypes"),
		new VCFInfoHeaderLine(MINOR_ALLELE_FREQUENCY_KEY, 1, VCFHeaderLineType.Float, "Minor Allele Frequency based in individual sample genotypes")));
	
	

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			AnnotatorCompatibleWalker walker, ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		
		int[] alleleFrequencies;
		int totalFrequency;
		
		char[] ac = initializeAlleleCharsArray();
		int[] ai = initializeAlleleIndecesArray();
		int alleleCount = initializeAlleleArrays(vc, ac, ai);
		
		alleleFrequencies = new int[alleleCount];
        totalFrequency = 0;
		for (Genotype genotype : vc.getGenotypes() ) {
			int plodicity = genotype.getPloidy();
			for (int i = 0; i < plodicity; i++) {
				Allele a = genotype.getAllele(i);
				if (a.isNoCall())
					continue;
				int aindex = ai[a.getBases()[0] - Byte.MIN_VALUE];
				alleleFrequencies[aindex]++;
			}
			totalFrequency += plodicity;
		}
		int minIndex = 0;
		int minFrequency = alleleFrequencies[ai[ac[0] - Byte.MIN_VALUE]];
		for (int i = 1; i < alleleCount; i++) {
			int aindex = ai[ac[i] - Byte.MIN_VALUE];
			int afreq = alleleFrequencies[aindex];
			if (afreq < minFrequency) minIndex = i;
		}

        Map<String, Object> map = new HashMap<String,Object>(2);
        map.put(MINOR_ALLELE_KEY, ac[minIndex] );
        map.put(MINOR_ALLELE_FREQUENCY_KEY, formatFractionAnnotation(((double) minFrequency) / (double) totalFrequency));
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
	
	protected static final int NO_ALLELE = -1;
	
	protected ThreadLocal<int[]> alleleIndeces = new ThreadLocal<int[]>();
	protected ThreadLocal<char[]> alleleChars = new ThreadLocal<char[]>();

	protected int initializeAlleleArrays(VariantContext vc, char[] ac, int[] ai) {
		int nextIndex = 0;
		Allele refAllele = vc.getReference();
		ac[nextIndex] = refAllele.getBaseString().charAt(0);
		ai[vc.getReference().getBases()[0] - Byte.MIN_VALUE] = nextIndex++;
		for (Allele a : vc.getAlternateAlleles()) {
			ac[nextIndex] = a.getBaseString().charAt(0);
			ai[a.getBases()[0] - Byte.MIN_VALUE] = nextIndex++;
		}
		return nextIndex;
	}

	protected int[] initializeAlleleIndecesArray() {
		int[] ai = alleleIndeces.get();
		if (ai == null)
			alleleIndeces
					.set(ai = new int[Byte.MAX_VALUE - Byte.MIN_VALUE + 1]);
		else
			Arrays.fill(ai, NO_ALLELE);
		return ai;
	}

	protected char[] initializeAlleleCharsArray() {
		char[] ac = this.alleleChars.get();
		if (ac == null)
			this.alleleChars.set(ac = new char[Byte.MAX_VALUE - Byte.MIN_VALUE
					+ 1]);
		return ac;
	}


}
