package net.malariagen.gatk.walker;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import net.malariagen.utils.Nucleotide;
import net.malariagen.utils.Trinucleotide;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;

public class SequenceComplexity {

	private static final boolean doTriEnt = false;
	int windowSize;
	int size;
	List<GenomeLoc> locs;
	GenomeLoc end;
	LinkedList<Nucleotide> nucs;
	LinkedList<Trinucleotide> trinucs;
	LinkedList<Integer> refMQs;
	int[] nucsCount;
	int[] trinucsCount;
	int nucsTotal;
	int trinucsTotal;
	private double[] prob;
	private double[] inverseProbLog;
	private double[] entropy;
	private double[] entropyMinus2;

	public List<LocusComplexity> flush() {
		List<LocusComplexity> result = new ArrayList<LocusComplexity>(nucs.size());
		while (nucs.size() > 1) {
		  removeNucleotide();
                  result.add(emit());
		}
                if (nucs.size() > 0) removeNucleotide();
		return result;
	}

	
	public LocusComplexity count(ReferenceContext ref,int refMQ) {
		return count(ref,refMQ,false);
	}
	
	public LocusComplexity count(ReferenceContext ref,int refMQ, boolean forceEmit) {
		GenomeLoc loc = ref.getLocus();
		Nucleotide n = Nucleotide.fromByte(ref.getBase());
		if (end != null && !loc.getContig().equals(end.getContig()))
			clear();
		int gap = (nucs.size() != 0) ? end.getStop() - loc.getStart() : 0;
		while (gap-- > 0) {
			if (nucs.size() == windowSize) 
				removeNucleotide();
			trinucs.add(Trinucleotide.NNN);
			refMQs.add(refMQ);
			nucs.add(Nucleotide.N);
			locs.add(loc);
			loc = ref.getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart() + 1, loc.getStart() + 1);
		}
		if (nucs.size() == windowSize) 
			removeNucleotide();
		end = loc;
		nucs.add(n);
		refMQs.add(refMQ);
		locs.add(loc);
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
		trinucs.add(t);

		if (t != Trinucleotide.NNN) {
			trinucsCount[t.ordinal()]++;
			trinucsTotal++;
		}

		if (forceEmit || nucs.size() == windowSize)
			return emit();
		else if (nucs.size() > windowSize)
			throw new RuntimeException("cannot be " + windowSize + "  "  + nucs.size());
		else
			return null;

	}


	private void removeNucleotide() {
		Nucleotide n0 = nucs.remove(0);
		Trinucleotide t0 = trinucs.remove(0);
		refMQs.remove(0);
		locs.remove(0);

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
	}

	// Fix seed to make it deterministic.
	private static Random rnd = new Random(13);

	private LocusComplexity emit() {
		GenomeLoc loc = locs.get(0);
		double gcBias = Double.NaN;
		double nucEnt = 0;
		double gcHet = Double.NaN;
		double triEnt = 0;
		int[] nucsCount = this.nucsCount;
		int[] trinucsCount = this.trinucsCount;
		// if there is some ambiguous nucs in the window we adjust the counts up as to emulate the missing one
		// keeping the proportion seen in the others.
		if (windowSize != nucsTotal) {
			int delta = windowSize - nucsTotal;
			int triDelta = windowSize - 2 - trinucsTotal;
			nucsCount = Arrays.copyOf(nucsCount, nucsCount.length);
			trinucsCount = Arrays.copyOf(trinucsCount, trinucsCount.length);
			double invNucsTotal = 1.0 / (double) nucsTotal;
			for (int i = 0; i < delta; i++) {
				double draw = rnd.nextDouble();
				for (int j = 0; j < 4; j++)
					if (j == 3
							|| (draw -= invNucsTotal * this.nucsCount[j]) <= 0) {
						nucsCount[j]++;
						break;
					}
				if (i >= triDelta) continue; 
				draw = rnd.nextDouble();
				if (doTriEnt)
				for (int j = 0; j < 64; j++)
					if (j == 63
							|| (draw -= invNucsTotal * this.trinucsCount[j]) <= 0) {
						trinucsCount[j]++;
						break;
					}
			}
		}

		// if (windowSize == nucsTotal) {
		gcBias = (nucsCount[1] + nucsCount[2]) * prob[1];
		gcBias *= 100;
		for (int i = 0; i < 4; i++)
			nucEnt += entropy[nucsCount[i]];
		if (doTriEnt) 
			for (int i = 0; i < 64; i++)
			  triEnt += entropyMinus2[trinucsCount[i]];
		LocusComplexity result = new LocusComplexity(loc, nucsTotal, nucs.get(0), gcBias, nucEnt,
				gcHet, triEnt);
		result.setRefMQ(refMQs.get(0));
		return result;
	}

	SequenceComplexity(int ws) {
		if (ws <= 0)
			throw new IllegalArgumentException(
					"invalid window size must be greater than 0");
		this.windowSize = ws;
		prob = new double[ws + 1];
		inverseProbLog = new double[ws + 1];
		entropy = new double[ws + 1];
		entropyMinus2 = new double[ws - 1];
		nucs = new LinkedList<Nucleotide>();
		refMQs = new LinkedList<Integer>();
		trinucs = new LinkedList<Trinucleotide>();
		nucsCount = new int[4];
		trinucsCount = new int[64];
		locs = new LinkedList<GenomeLoc>();
		prob[0] = 0.5 / (double) ws;
		inverseProbLog[0] = Math.log(prob[0]);
		entropy[0] = 0;
		entropyMinus2[0] = 0;
		double logWsMinus2 = Math.log(ws - 2);

		for (int i = 1; i <= ws; i++) {
			prob[i] = i / (double) ws;
			inverseProbLog[i] = Math.log(prob[i]);
			entropy[i] = -prob[i] * inverseProbLog[i];
			if (i <= (ws - 2))
				entropyMinus2[i] = -(((double) i) / (double) (ws - 2))
						* (Math.log(i) - logWsMinus2);
		}
	}

	public static SequenceComplexity create(int windowSize) {
		return new SequenceComplexity(windowSize);
	}

	public class LocusComplexity {

		private double gcBias;
		private double nucEnt;
		private double gcHet;
		private GenomeLoc loc;
		private Nucleotide refNuc;
		private double triEnt;
		private int size;
		private int refMQ;

		public LocusComplexity(GenomeLoc loc, int size, Nucleotide n,
				double gcBias, double nucEnt, double gcHet, double triEnt) {
			this.loc = loc;
			this.gcBias = gcBias;
			this.nucEnt = nucEnt;
			this.gcHet = gcHet;
			this.refNuc = n;
			this.triEnt = triEnt;
			this.size = size;
		}

		public int size() {
			return size;
		}

		public GenomeLoc getLocus() {
			return loc;
		}

		public Nucleotide getRefNuc() {
			return refNuc;
		}

		public double getGcBias() {
			return gcBias;
		}

		public double getNucEnt() {
			return this.nucEnt;
		}

		public double getGcHet() {
			return this.gcHet;
		}

		public double getTriEnt() {
			return this.triEnt;
		}

		public int getRefMQ() {
			return this.refMQ;
		}
		
		public void setRefMQ(int value) {
			this.refMQ = value;
		}

		public int getNucCount() {
			return size;
		}

	}


}
