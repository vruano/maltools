package net.malariagen.gatk.walker;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import net.malariagen.utils.Nucleotide;
import net.malariagen.utils.Trinucleotide;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;

public class SequenceComplexity {

	private static final int NUC_C_ORDINAL = Nucleotide.C.ordinal();
	private static final int NUC_G_ORDINAL = Nucleotide.G.ordinal();
	private int windowSize;
	private int windowSizeMinus2;
	int size;
	List<GenomeLoc> locs;
	GenomeLoc end;
	private LinkedList<Nucleotide> nucs;
	private LinkedList<Boolean> isCoding;
	LinkedList<Trinucleotide> trinucs;
	int[] nucsCount;
	int[] trinucsCount;
	int codingCount;
	int nucsTotal;
	int trinucsTotal;
	private double[] entropy;
	private double[] entropyMinus2;
	public int triSize;

	public List<LocusComplexity> flush() {
		List<LocusComplexity> result = new ArrayList<LocusComplexity>(nucs.size());
		while (locs.size() > 1) {
		   removeNucleotide();
//		   if (windowSize == 256) System.err.println("" + locs.get(0) + " " + windowSize);
		   result.add(emit());
		}
		clear();
		return result;
	}

	
	public LocusComplexity count(ReferenceContext ref) {
		return count(ref,null);
	}
	
	public LocusComplexity count(ReferenceContext ref, Boolean isCoding) {
		GenomeLoc loc = ref.getLocus();
		Nucleotide n = Nucleotide.fromByte(ref.getBase());
		if (end != null && !loc.getContig().equals(end.getContig()))
			clear();
		int gap = (end != null) ? end.getStop() - loc.getStart() : 0;
		while (gap-- > 0) {
			if (nucs.size() == windowSize) 
				removeNucleotide();
			trinucs.add(Trinucleotide.NNN);
			nucs.add(Nucleotide.N);
			this.isCoding.add(null);
			locs.add(loc);
			loc = ref.getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart() + 1, loc.getStart() + 1);
		}
		if (nucs.size() == windowSize) 
			removeNucleotide();
		end = loc;
		nucs.add(n);
		locs.add(loc);
		this.isCoding.add(isCoding);
		if (isCoding == Boolean.TRUE) codingCount++;
		if (n != Nucleotide.N) {
			nucsCount[n.ordinal()]++;
			nucsTotal++;
		} 

		Trinucleotide t = Trinucleotide.NNN;
		if (trinucs.size() > 0) {
			Trinucleotide last = trinucs.getLast();
			t = last.shift(n);
		}
		if (t == Trinucleotide.NNN && nucs.size() >= 3) {
			t = Trinucleotide.fromNucleotides(nucs.get(nucs.size() - 3),
					nucs.get(nucs.size() - 2), n);
		}
		
		if (nucs.size() > 2)
		  trinucs.add(t);

		if (t != Trinucleotide.NNN) {
			if (++trinucsCount[t.ordinal()] > windowSizeMinus2)
				throw new IllegalStateException("cannot be " + trinucs.size() + " " + nucs.size() + " " + windowSize + " " + trinucsCount[t.ordinal()] + " " + t + " " + nucsTotal + " " + trinucsTotal + " "  + Arrays.toString(nucs.toArray()) + " " + Arrays.toString(trinucs.toArray()));
			trinucsTotal++;
		}

		if (nucs.size() == windowSize)
			return emit();
		else if (nucs.size() > windowSize)
			throw new RuntimeException("cannot be " + windowSize + "  "  + nucs.size());
		else
			return null;

	}


	private void removeNucleotide() {

		Nucleotide n0 = nucs.remove(0);
		Trinucleotide t0 = !trinucs.isEmpty() ? trinucs.remove(0) : Trinucleotide.NNN;
		locs.remove(0);
		Boolean coding = isCoding.remove();
		if (coding == Boolean.TRUE)
			codingCount--;
		if (n0 != Nucleotide.N) {
			nucsCount[n0.ordinal()]--;
			nucsTotal--;
		}
		if (t0 != Trinucleotide.NNN) {
			trinucsCount[t0.ordinal()]--;
			trinucsTotal--;
		}
	}

	private void clear() {
		end = null;
		nucs.clear();
		trinucs.clear();
		locs.clear();
		nucsTotal = 0;
		trinucsTotal = 0;
		Arrays.fill(nucsCount, 0);
		Arrays.fill(trinucsCount, 0);
		codingCount = 0;
		isCoding.clear();
	}

	// Fix seed to make it deterministic.
//	private static Random rnd = new Random(13);

	private LocusComplexity emit() {
		GenomeLoc loc = locs.get(0);
		int gcCount = 0;
		double nucEnt = 0;
		double triEnt = 0;
		double dust = 0;
		gcCount = (nucsCount[NUC_G_ORDINAL] + nucsCount[NUC_C_ORDINAL]);
		nucEnt = entropy(nucsCount,nucsTotal);
		triEnt = entropy(trinucsCount,trinucsTotal);
		for (int c : trinucsCount) 
			if (c > 0) dust += (c * (c-1)) >> 1;
		if (trinucsTotal < windowSizeMinus2)
			dust = (dust * (windowSizeMinus2)) / ((double) trinucsTotal * (trinucsTotal - 1));
		else 
			dust = dust / (trinucsTotal - 1);
		
		LocusComplexity result = new LocusComplexity(loc, nucsTotal, trinucsTotal, nucs.get(nucs.size() - locs.size()), gcCount, nucEnt,
				dust, triEnt, codingCount);
		return result;
	}

	private double entropy(int[] counts, int total) {
		double result = 0;
		if (total == windowSize) 
			for (int c : counts)
				result += entropy[c];
		else if (total == windowSizeMinus2) 
			for (int c : counts)
				result += entropyMinus2[c];
		else 
			for (int c : counts) {
				if (c == 0) continue;
				double prob = ((double) c) / ((double)total); 
				result -= prob * Math.log(prob);
			}
		return result;
	}


	SequenceComplexity(int ws) {
		if (ws <= 0)
			throw new IllegalArgumentException(
					"invalid window size must be greater than 0");
		windowSize = ws;
		windowSizeMinus2 = ws - 2;
		double[] prob = new double[ws + 1];
		double[] lnProb = new double[ws + 1];
		entropy = new double[ws + 1];
		entropyMinus2 = new double[ws - 1];
		nucs = new LinkedList<Nucleotide>();
		isCoding = new LinkedList<Boolean>();
		trinucs = new LinkedList<Trinucleotide>();
		nucsCount = new int[4];
		trinucsCount = new int[64];
		locs = new LinkedList<GenomeLoc>();
		prob[0] = 0;
		lnProb[0] = Math.log(prob[0]);
		entropy[0] = 0;
		entropyMinus2[0] = 0;
		double logWsMinus2 = Math.log(windowSizeMinus2);

		for (int i = 1; i <= ws; i++) {
			prob[i] = i / (double) ws;
			lnProb[i] = Math.log(prob[i]);
			entropy[i] = -prob[i] * lnProb[i];
			if (i <= windowSizeMinus2)
				entropyMinus2[i] = -(((double) i) / (double) windowSizeMinus2)
						* (Math.log(i) - logWsMinus2);
		}
	}

	public static SequenceComplexity create(int windowSize) {
		return new SequenceComplexity(windowSize);
	}
	
	public static SequenceComplexity create(int windowSize, int padding) {
		SequenceComplexity result = create(windowSize);
		result.pad(padding);
		return result;
	}

	private void pad(int padding) {
		for (int i = 0; i < padding; i++) {
			nucs.add(Nucleotide.N);
			isCoding.add(null);
			if (i > 1)
			  trinucs.add(Trinucleotide.NNN);
		}
	}

	public class LocusComplexity {

		private double gcCount;
		private double nucEnt;
		private double dust;
		private GenomeLoc loc;
		private Nucleotide refNuc;
		private double triEnt;
		private int size;
		private int triSize;
		private int codingCount;

		public LocusComplexity(GenomeLoc loc, int size, int triSize, Nucleotide n,
				int gcCount, double nucEnt, double dust, double triEnt, int codingCount) {
			this.loc = loc;
			this.gcCount = gcCount;
			this.nucEnt = nucEnt;
			this.dust = dust;
			this.refNuc = n;
			this.triEnt = triEnt;
			this.size = size;
			this.triSize = triSize;
			this.codingCount = codingCount;
		}

		public int size() {
			return size;
		}
		
		public int getCodingCount() {
			return codingCount;
		}

		public GenomeLoc getLocus() {
			return loc;
		}

		public Nucleotide getRefNuc() {
			return refNuc;
		}

		public double getGcCount() {
			return gcCount;
		}

		public double getNucEnt() {
			return this.nucEnt;
		}

		public double getDUST() {
			return this.dust;
		}

		public double getTriEnt() {
			return this.triEnt;
		}

		public int getNucCount() {
			return size;
		}
		
		public int getTrinucCount() {
			return triSize;
		}

	}


}
