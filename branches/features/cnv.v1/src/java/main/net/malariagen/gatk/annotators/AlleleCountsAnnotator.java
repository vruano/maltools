package net.malariagen.gatk.annotators;

import java.util.Arrays;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;



// Contains common methods that apply to annotators that need to do integer counting upon alleles.
public abstract class AlleleCountsAnnotator {

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
