package net.malariagen.utils;

import java.util.Arrays;

/**
 * Enumeration IUPAC code for nucleotides. 
 * 
 * @author Valentin Ruano-Rubio &lt;valentin.ruano@gmail.com&gt;
 *
 */
public enum NucleotideIUPAC {

	/**
	 * Concrete code to represent adenine.
	 */
	A(Mask.A),

	/**
	 * Concrete code to represent cytosine.
	 */
	C(Mask.C),

	/**
	 * Concrete code to represent guanine.
	 */
	G(Mask.G),

	/**
	 * Concrete code to represent thymine.
	 * <p/>
	 * Notice that this code is considered semantically equivalent to {@link #U}
	 * (see {@link #isEqualTo(NucleotideIUPAC) isEqualTo}).
	 */
	T(Mask.T),

	/**
	 * Concrete code to represent uracil bases: <i>U</i>.
	 * <p/>
	 * Notice that this code is considered semantically equivalent to {@link #T}
	 * (see {@link #isEqualTo(NucleotideIUPAC) isEqualTo}).
	 */
	U(Mask.U),

	/**
	 * Ambiguous code to represent purine bases: <i>A</i> or <i>G</i>
	 */
	R(Mask.A | Mask.G),

	/**
	 * Ambiguous code to represent pyrimidine bases: <i>C</i> or <i>T/U</i>
	 */
	Y(Mask.C | Mask.T),

	/**
	 * Ambiguous code to represent bases <i>G</i> or <i>C</i>
	 */
	S(Mask.G | Mask.C),

	/**
	 * Ambiguous code to represent bases <i>A</i> or <i>T/U</i>
	 */
	W(Mask.A | Mask.T),

	/**
	 * Ambiguous code to represent bases <i>A</i> or <i>C</i>
	 */
	K(Mask.G | Mask.T),

	/**
	 * Ambiguous code to represent bases <i>A</i> or <i>C</i>
	 */
	M(Mask.A | Mask.C),

	/**
	 * Ambiguous code to represent bases <i>C, G</i> or <i>T/U</i>.
	 */
	B(Mask.C | Mask.G | Mask.T),

	/**
	 * Ambiguous code to represent bases <i>A, G</i> or <i>T/U</i>.
	 */
	D(Mask.A | Mask.G | Mask.T),

	/**
	 * Ambiguous code to represent bases <i>A, C</i> or <i>T/U</i>.
	 */
	H(Mask.A | Mask.C | Mask.T),

	/**
	 * Ambiguous code to represent bases <i>A, C</i> or <i>G</i>.
	 */
	V(Mask.A | Mask.C | Mask.G),

	/**
	 * Ambiguous code to represent any base: <i>A, C, G</i> or <i>T/U</i>.
	 */
	N(Mask.A | Mask.C | Mask.G | Mask.T),

	/**
	 * Special code that represents a gap. Although '.' and '-' are both valid
	 * IUPAC representations of a gap, '-' takes preference in this enumeration
	 * and is the one used in values returned by this implementation
	 */
	GAP(16, (byte) '-'),

	/**
	 * Special code to represent, in fact, the lack of one. NaC stands for
	 * <i>"Not a Code"</i>.
	 */
	NaC(32);

	class Mask {
		public static final int A = 1;
		public static final int C = 2;
		public static final int G = 4;
		public static final int T = 8;
		public static final int U = 8;
	}

	private byte byteValue = -1;

	private int maskValue;

	private static final NucleotideIUPAC[] BYTE_TO_CODE = new NucleotideIUPAC[128];

	private static final NucleotideIUPAC[] MASK_TO_CODE = new NucleotideIUPAC[17];

	static {
		Arrays.fill(BYTE_TO_CODE, NaC);
		Arrays.fill(MASK_TO_CODE, NaC);
		for (NucleotideIUPAC i : NucleotideIUPAC.values()) {
			byte lowerCase = (byte) Character.toLowerCase(i.byteValue);
			if (i.byteValue >= 0)
				BYTE_TO_CODE[i.byteValue] = BYTE_TO_CODE[lowerCase] = i;
			if (i.maskValue < 16)
				MASK_TO_CODE[i.maskValue] = i;
		}
		BYTE_TO_CODE[(byte) '.'] = GAP;
	}

	NucleotideIUPAC(int mask) {
		maskValue = mask;
		byteValue = (byte) this.toString().charAt(0);
	}

	NucleotideIUPAC(int mask, byte b) {
		maskValue = mask;
		byteValue = b;
	}

	/**
	 * 
	 * @param b
	 *            that represent the requested code.
	 * @return Never {@code null} but {@link #NaC} if byte is not a legal IUPAC
	 *         code byte.
	 */
	public static NucleotideIUPAC fromBase(byte b) {
		if (b < 0)
			return NaC;
		return BYTE_TO_CODE[b];
	}

	/**
	 * Checks whether this code is compatible with a given base.
	 * 
	 * @param b
	 *            the byte representation of the base.
	 * @return {@code true} if this code can represent such a base,
	 *         {@code false} otherwise.
	 */
	public boolean isCompatibleWith(byte b) {
		NucleotideIUPAC iupac = fromBase(b);
		return iupac == null ? false : isCompatibleWith(iupac);
	}

	/**
	 * Checks whether this and another code has some possible base in common in
	 * which case they are considered compatible. Also {@link #GAP} is
	 * compatible with and only with itself.
	 * 
	 * <p/>
	 * 
	 * Notice however that {@link #NaC} is incompatible with any code including
	 * itself.
	 * 
	 * @param o
	 *            the other code to compare with.
	 * @throws NullPointerException
	 *             if {@code o} is {@code null}.
	 * @return {@code true} if these are compatible, {@code false} otherwise.
	 */
	public boolean isCompatibleWith(NucleotideIUPAC o) {
		if (o == NaC || this == NaC)
			return false;
		else if (this == o)
			return true;
		else
			return (o.maskValue & this.maskValue) != 0;
	}

	/**
	 * Checks whether to (ambiguous) bases are compatible with each other.
	 * 
	 * <p/>
	 * This is equivalent to:
	 * <p/>
	 * {@code fromBase(b1).isCompatibleWith(fromBase(b2)) }
	 * 
	 * @param b1
	 *            first base to compare
	 * @param b2
	 *            second base to compare
	 * @return {@code true} if both are mutually compatible, {@code false}
	 *         otherwise.
	 */
	public static boolean areCompatible(byte b1, byte b2) {
		NucleotideIUPAC i1 = fromBase(b1);
		return i1 == null ? false : i1.isCompatibleWith(b2);
	}

	/**
	 * Returns the least ambiguous code that is compatible with all bases
	 * provided as an input.
	 * 
	 * @param b
	 *            set of bases to encode.
	 * @return never {@code null} but perhaps {@link #NaC} if {@code b} has zero
	 *         lenght or contains no-sense bases or a combination of nucleotide
	 *         and gap code bytes.
	 */
	public static NucleotideIUPAC fromBases(byte... b) {
		if (b.length == 0)
			return NaC;
		NucleotideIUPAC result = fromBase(b[0]);
		for (int i = 1; i < b.length; i++)
			result = result.union(fromBase(b[i]));
		return result;
	}

	/**
	 * Returns the code that represents the union between this and the other
	 * code provided.
	 * <p/>
	 * For sense no-gap code combination this would result in the most specific
	 * code that is compatible with all possible nucleotides based in the input
	 * codes.
	 * 
	 * @param o
	 *            the other code to combine with
	 * @return never {@code null} but {@link #NaC} when this or {@code o} is
	 *         {@link #NaC} or they result in a combination of gap and base codes.
	 */
	private NucleotideIUPAC union(NucleotideIUPAC o) {
		int unionMask = this.maskValue | o.maskValue;
		if (unionMask >= MASK_TO_CODE.length)
			return NaC;
		else
			return MASK_TO_CODE[unionMask];
	}

	public static boolean areEqual(byte b1, byte b2) {
		return fromBase(b1).isEqualTo(fromBase(b2));
	}

	/**
	 * Semantics correct comparison for two codes.
	 * 
	 * <p/>
	 * 
	 * In most case you will want to use this method instead of the ordinary
	 * {@link #equals(Object) equals} that is bound to strictly follow the Java
	 * enumeration semantics.
	 * 
	 * <p/>
	 * 
	 * This method and {@link #equals(Object) equals} are mostly equivalent
	 * except for:
	 * <ul>
	 * <li> {@link #T} is equal to {@link #U} and vice versa,</li>
	 * <li> {@link #NaC} is not equal to everything including itself (thus
	 * {@code NaC.isEqualTo(NaC) == false}).</li>
	 * </ul>
	 * 
	 * @param o
	 *            the other code to compare with.
	 * @throws NullPointerException
	 *             if {@code o} is {@code null}
	 * @return {@code true} if this and {@code o} are equal to each other,
	 *         {@code false} otherwise
	 */
	public boolean isEqualTo(NucleotideIUPAC o) {
		return false;
	}

	/**
	 * Returns the byte representation of this code.
	 * <p/>
	 * It case of letters it will return the uppercase form. For gaps it returns
	 * '-' over the alternative form '.'. For {@link #NaC} this is {@code -1}.
	 */
	public byte byteValue() {
		return byteValue;
	}

}
