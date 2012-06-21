package net.malariagen.gatk.coverage;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.malariagen.utils.Nucleotide;
import net.malariagen.utils.Trinucleotide;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;

public class SequenceComplexity {

	int windowSize;
	int size;
	List<GenomeLoc> locs;
	GenomeLoc end;
	LinkedList<Nucleotide> nucs;
	LinkedList<Trinucleotide> trinucs;
	int[] nucsCount;
	int[] trinucsCount;
	int nucsTotal;
	private double[] prob;
	private double[] inverseProbLog;
	private double[] entropy;
	private double[] entropyMinus2;

	public LocusComplexity count(ReferenceContext ref) {
		GenomeLoc loc = ref.getLocus();
		Nucleotide n = Nucleotide.fromByte(ref.getBase());
		int gap = (nucs.size() != 0) ? end.getStop() - loc.getStart() : 0;
		while (gap++ > 0) {
			if (nucs.size() == windowSize) {
				Nucleotide n0 = nucs.remove(0);
				Trinucleotide t0 = trinucs.remove(0);
				locs.remove(0);
				if (n0 != Nucleotide.NaN) {
					nucsCount[n0.ordinal()]--;
					nucsTotal--;
				}
				if (t0 != Trinucleotide.NaT)
					trinucsCount[t0.ordinal()]--;
			}
			trinucs.add(Trinucleotide.NaT);
			nucs.add(Nucleotide.NaN);
			locs.add(null);
		}
		if (nucs.size() == windowSize) {
			Nucleotide n0 = nucs.remove(0);
			Trinucleotide t0 = trinucs.remove(0);
			locs.remove(0);
			
			if (n0 != Nucleotide.NaN) {
				nucsCount[n0.ordinal()]--;
				nucsTotal--;
			}
			if (t0 != Trinucleotide.NaT)
				trinucsCount[t0.ordinal()]--;
			
		}
		end = loc;
		nucs.add(n);
		if (n != Nucleotide.NaN) {
			locs.add(loc);
			nucsCount[n.ordinal()]++;
			nucsTotal++;
		}
		else 
			locs.add(null);

		Trinucleotide t = Trinucleotide.NaT;
		if (trinucs.size() > 0) {
			Trinucleotide last = trinucs.getLast();
			t = last.shift(n);
		}
		if (t == Trinucleotide.NaT && nucs.size() >= 3) {
			t = Trinucleotide.fromNucleotides(nucs.get(nucs.size() - 3),nucs.get(nucs.size() - 2), n);
		}
		trinucs.add(t);
		if (t != Trinucleotide.NaT) 
			trinucsCount[t.ordinal()]++;
	
		
		if (nucsTotal == windowSize)
			return emit();
		else
			return null;

	}

	private LocusComplexity emit() {
//		if (true) throw new RuntimeException("Emit called");
		GenomeLoc loc = locs.get(0);
		double gcBias = Double.NaN;
		double nucEnt = 0;
		double gcHet = Double.NaN;

		gcBias = (nucsCount[1] + nucsCount[2]) * prob[1];
		for (int i = 0; i < 4; i++)
			nucEnt += entropy[nucsCount[i]];

		double maxHet = 0.5 * ( windowSize * windowSize - windowSize);

		int het = 0;
		int gc = 0;
		Iterator<Nucleotide> it = nucs.iterator();
		for (int i = 0; i < windowSize; i++) {
			Nucleotide n = it.next();
			switch (n) {
			case C:
			case G:
				gc++;
			default:
			}
			double expected = gcBias * (i + 1);
			if (expected > gc)
				het += Math.floor(expected) - gc;
			else
				het += gc - Math.ceil(expected);
		}
		gcHet = ((double)het) / (double) maxHet;

		gcBias *= 100;

		return new LocusComplexity(loc, nucs.get(0), gcBias, nucEnt, gcHet);
	}

	SequenceComplexity(int ws) {
		if (ws <= 0)
			throw new IllegalArgumentException(
					"invalid window size must be greater than 0");
		this.windowSize = ws;
		prob = new double[ws];
		inverseProbLog = new double[ws];
		entropy = new double[ws];
		entropyMinus2 = new double[ws - 2];
		nucs = new LinkedList<Nucleotide>();
		nucsCount = new int[4];
		trinucsCount = new int[64];
		locs = new LinkedList<GenomeLoc>();
		prob[0] = 0.5 / (double) ws;
		inverseProbLog[0] = Math.log(prob[0]);
		entropy[0] = 0;
		entropyMinus2[0] = 0;
		double logWsMinus2 = Math.log(ws - 2);
		
		for (int i = 1; i < ws; i++) {
			prob[i] = i / (double) ws;
			inverseProbLog[i] = Math.log(prob[i]);
			entropy[i] = -prob[i] * inverseProbLog[i];
			if (i < (ws - 2)) 
				entropyMinus2[i] = i / (double) (ws - 2)  * (Math.log(i) - logWsMinus2); 
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

		public LocusComplexity(GenomeLoc loc, Nucleotide n, double gcBias,
				double nucEnt, double gcHet) {
			this.loc = loc;
			this.gcBias = gcBias;
			this.nucEnt = nucEnt;
			this.gcHet = gcHet;
			this.refNuc = n;
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

	}

}
