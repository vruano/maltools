package net.malariagen.gatk.coverage;


import java.util.Arrays;
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
	int trinucsTotal;
	private double[] prob;
	private double[] inverseProbLog;
	private double[] entropy;
	private double[] entropyMinus2;

	public LocusComplexity count(ReferenceContext ref) {
		GenomeLoc loc = ref.getLocus();
		Nucleotide n = Nucleotide.fromByte(ref.getBase());
		if (end != null && !loc.getContig().equals(end.getContig()))
			clear();
		int gap = (nucs.size() != 0) ? end.getStop() - loc.getStart() : 0;
		while (gap++ > 0) {
			if (nucs.size() == windowSize) {
				Nucleotide n0 = nucs.remove(0);
				Trinucleotide t0 = trinucs.remove(0);
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
			trinucs.add(Trinucleotide.NNN);
			nucs.add(Nucleotide.N);
			locs.add(loc);
		}
		if (nucs.size() == windowSize) {
			Nucleotide n0 = nucs.remove(0);
			Trinucleotide t0 = trinucs.remove(0);
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
		end = loc;
		nucs.add(n);
		if (n != Nucleotide.N) {
			locs.add(loc);
			nucsCount[n.ordinal()]++;
			nucsTotal++;
		}
		else 
			locs.add(loc);

		Trinucleotide t = Trinucleotide.NNN;
		if (trinucs.size() > 0) {
			Trinucleotide last = trinucs.getLast();
			t = last.shift(n);
		}
		if (t == Trinucleotide.NNN && nucs.size() >= 3) {
			t = Trinucleotide.fromNucleotides(nucs.get(nucs.size() - 3),nucs.get(nucs.size() - 2), n);
		}
		trinucs.add(t);
		if (t != Trinucleotide.NNN) { 
			trinucsCount[t.ordinal()]++;
			trinucsTotal++;
		}
	
		
		if (nucs.size() == windowSize)
			return emit();
		else
			return null;

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

	private LocusComplexity emit() {
//		if (true) throw new RuntimeException("Emit called");
		GenomeLoc loc = locs.get(0);
		double gcBias = Double.NaN;
		double nucEnt = 0;
		double gcHet = Double.NaN;
		double triEnt = 0;

		if (windowSize == nucsTotal) { 
		  gcBias = (nucsCount[1] + nucsCount[2]) * prob[1];
		  gcBias *= 100;
		  for (int i = 0; i < 4; i++)
			  nucEnt += entropy[nucsCount[i]];
		  for (int i = 0; i < 64; i++)
			  triEnt += entropyMinus2[trinucsCount[i]];
		}	
		else {
		  double invNucsTotal =  1.0 / (double) nucsTotal;
		  double logNucsTotal = Math.log(nucsTotal);
		  gcBias = (nucsCount[1] + nucsCount[2]) * invNucsTotal;
		  gcBias *= 100;
		  for (int i = 0; i < 4; i++)
			  nucEnt -= nucsCount[i] == 0 ? 0 : nucsCount[i] * invNucsTotal * (Math.log(nucsCount[i]) - logNucsTotal);
		  double invTrinucsTotal = 1.0 / (double) trinucsTotal;
		  double logTrinucsTotal = Math.log(trinucsTotal);
		  for (int i = 0; i < 64; i++)
			  triEnt -= trinucsCount[i] == 0 ? 0 : trinucsCount[i] * invTrinucsTotal * (Math.log(trinucsCount[i]) - logTrinucsTotal);		  
		}
		return new LocusComplexity(loc, nucsTotal, nucs.get(0), gcBias, nucEnt, gcHet, triEnt);
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
		trinucs = new LinkedList<Trinucleotide>();
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
				entropyMinus2[i] = - (((double) i) / (double) (ws - 2))  * (Math.log(i) - logWsMinus2); 
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

		public LocusComplexity(GenomeLoc loc, int size, Nucleotide n, double gcBias,
				double nucEnt, double gcHet, double triEnt) {
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

	}

}
