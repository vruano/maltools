package net.malariagen.utils;

public enum Trinucleotide {
	AAA, AAC, AAG, AAT,
	ACA, ACC, ACG, ACT,
	AGA, AGC, AGG, AGT,
	ATA, ATC, ATG, ATT,
	CAA, CAC, CAG, CAT,
	CCA, CCC, CCG, CCT,
	CGA, CGC, CGG, CGT,
	CTA, CTC, CTG, CTT,
	GAA, GAC, GAG, GAT,
	GCA, GCC, GCG, GCT,
	GGA, GGC, GGG, GGT,
	GTA, GTC, GTG, GTT,
	TAA, TAC, TAG, TAT,
	TCA, TCC, TCG, TCT,
	TGA, TGC, TGG, TGT,
	TTA, TTC, TTG, TTT,
	NNN;
	
	
//	private static Trinucleotide[] BYTE_TO_TRI = new Trinucleotide[Byte.MAX_VALUE];
//
//	static {
//		Arrays.fill(BYTE_TO_TRI, NaT);
//		Trinucleotide[] values = values();
//		for (int i = 0; i < values.length; i++) 
//			BYTE_TO_TRI[i] = values[i];
//	}
//	
//	Trinucleotide fromByte(byte b) {
//		if (b < 0)
//			return NaT;
//		return BYTE_TO_TRI[b];
//	}
//	
	
	public Trinucleotide shift(Nucleotide n) {
		if (n == Nucleotide.N)
			return NNN;
		int o = ordinal();
		o = ((o << 2) + n.ordinal()) & 0x4F;
		return values()[o];
	}
	
	public Trinucleotide unshift(Nucleotide n) {
		if (n == Nucleotide.N)
			return NNN;
		int o = ordinal();
		o = (o >> 2) + (n.ordinal() << 4);
		return values()[o];
		
	}
	
	public static Trinucleotide fromNucleotides(Nucleotide n2, Nucleotide n1, Nucleotide n0) {
		if (n2 == Nucleotide.N || n1 == Nucleotide.N || n0 == Nucleotide.N)
			return NNN;
		int o = (n2.ordinal() << 4) + (n1.ordinal() << 2) + n0.ordinal();
		return values()[o];
	}
	
	public Nucleotide nucleotide(int p) {
		if (p < 0 || p > 2)
			throw new IllegalArgumentException("the postition to mutate must be 0, 1 or 2");
		int o = ordinal() >> ((2 - p) << 1);
		return Nucleotide.values()[o];
	}
	
	public Trinucleotide mutate(int p, Nucleotide n) {
		if (p < 0 || p > 2)
			throw new IllegalArgumentException("the postition to mutate must be 0, 1 or 2");
		if (n == Nucleotide.N)
			return NNN;
		int shift = (2 - p) << 1;
		int mask = 0x04 << shift;
		int o = (ordinal() & mask)  + n.ordinal() << shift;
		return Trinucleotide.values()[o];		
	}

}
