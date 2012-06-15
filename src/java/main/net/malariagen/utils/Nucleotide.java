package net.malariagen.utils;

import java.util.Arrays;

public enum Nucleotide {

	A, C, G, T, NoN;
	
	private static final Nucleotide[] BYTE_TO_CODE = new Nucleotide[Byte.MAX_VALUE];
	
	static {
		Arrays.fill(BYTE_TO_CODE, NoN);
		BYTE_TO_CODE[(byte)'a'] = BYTE_TO_CODE[(byte)'A'] = A;
		BYTE_TO_CODE[(byte)'c'] = BYTE_TO_CODE[(byte)'C'] = C;
		BYTE_TO_CODE[(byte)'g'] = BYTE_TO_CODE[(byte)'G'] = G;
		BYTE_TO_CODE[(byte)'t'] = BYTE_TO_CODE[(byte)'T'] = T;
		BYTE_TO_CODE[(byte)'u'] = BYTE_TO_CODE[(byte)'U'] = T;
	}
	
	public static Nucleotide fromByte(byte b) {
		if (b < 0)
			return NoN;
		return BYTE_TO_CODE[b];
	}
	
	public static Nucleotide fromIUPAC(NucleotideIUPAC iupac) {
		switch (iupac) {
		case A: return A;
		case C: return C;
		case G: return G;
		case T:
		case U: return T;
		default:
			return NoN;
		}
	}

	public NucleotideIUPAC toIUPAC() {
		switch (this) {
		case A: return NucleotideIUPAC.A;
		case C: return NucleotideIUPAC.C;
		case G: return NucleotideIUPAC.G;
		case T: return NucleotideIUPAC.T;
		default: 
			return NucleotideIUPAC.NaC;
		}
	}
	
	public byte byteValue() {
		switch (this) {
		case A: return (byte) 'A';
		case C: return (byte) 'C';
		case G: return (byte) 'G';
		case T: return (byte) 'T';
		default: return -1;
		}
	}
	
}
