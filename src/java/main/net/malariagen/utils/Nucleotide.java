package net.malariagen.utils;

import java.util.Arrays;

public enum Nucleotide {
	A, C, G, T, N;

	private final byte byteValue;
	
	private static final Nucleotide[] FROM_BYTE = new Nucleotide[Byte.MAX_VALUE];
	
	static {
		Arrays.fill(FROM_BYTE,N);
		FROM_BYTE['a'] = FROM_BYTE['A'] = Nucleotide.A;
		FROM_BYTE['c'] = FROM_BYTE['C'] = Nucleotide.C;
		FROM_BYTE['g'] = FROM_BYTE['G'] = Nucleotide.G;
		FROM_BYTE['t'] = FROM_BYTE['T'] = Nucleotide.T;
	}
	
	Nucleotide() {
		byteValue = (byte) toString().charAt(0);
	}
	
	public static Nucleotide fromByte(byte b) {
		if (b < 0)
			return N;
		return FROM_BYTE[b];
	}

	public byte byteValue() {
		return byteValue;
	}

	public boolean isProper() {
		return ordinal() <= T.ordinal();
	}
}
