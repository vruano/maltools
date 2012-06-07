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
	 * Checks whether this code is specific, that is it is compatible with a single nucleotide xor a gap.
	 * @return {@code true} if it is specific, {@code false} otherwise
	 */
	public boolean isUnambiguous() {
		switch(this) {
		case T:
		case A:
		case G:
		case C:
		case GAP:
			return true;
		default:
			return false;
		}
	}
	
	public static boolean  isUnambiguous(byte b) {
		return fromBase(b).isUnambiguous();
	}
	
	/**
	 * Checks whether this code is ambiougous, that is it is a nucleotide code and there is more than one nucleotide compatible.
	 * 
	 * @return {@code true} if it is ambigous, {@code false} otherwise.
	 */
	public boolean isAmbiguous() {
		switch (this) {
		case T:
		case A:
		case G:
		case C:
		case GAP:
		case NaC:
			return false;
		default:
			return true;
		}
	}

	public static boolean  isAmbiguous(byte b) {
		return fromBase(b).isAmbiguous();
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
	 * codes.
	 * <p/>
	 * For a nucleotide and and a non-nucleotide combination this results in
	 * {@link #NaC}. Thus for example,
	 * <code> {@link #N}.union({@link #GAP}) == {@link #NaC} </code>.
	 * <p/>
	 * Any union involving {@link #NaC} results in {@link #NaC}.
	 * <p/>
	 * In the special case that the union is {@link #T}, thus {@link #U}, the
	 * returned code, if at all different is the current one, so:
	 * <code> {@link #T}.union({@link #U}) == {@link #T}.union({@link #T}) == {@link #T} != </code>
	 * <code> {@link #U}.union({@link #U}) == {@link #U}.union({@link #T}) == {@link #U} </code>
	 * 
	 * @param o
	 *            the other code to combine with
	 * @return never {@code null} but {@link #NaC} when this or {@code o} is
	 *         {@link #NaC} or they result in a combination of gap and base
	 *         codes.
	 */
	public NucleotideIUPAC union(NucleotideIUPAC o) {
		int unionMask = this.maskValue | o.maskValue;
		if (unionMask == T.maskValue)
			return this;
		if (unionMask >= MASK_TO_CODE.length)
			return NaC;
		else
			return MASK_TO_CODE[unionMask];
	}

	/**
	 * In the special case that the intersection is equal to {@link #T}, thus
	 * {@link #U}, it returns one or the other strictly depending on the value
	 * given to {@code isDNA}. Notice that this is so even in the degenerated case where both codes are
	 * equal to either of these two. For example:
	 * <code> {@link #T}.intersection({@link #T},false) == {@link #U} </code>
	 */
	private NucleotideIUPAC intersection(NucleotideIUPAC o, boolean isDNA) {
	  throw new RuntimeException("to be implemented");
	}

	public static boolean areEqual(byte b1, byte b2) {
		return fromBase(b1).equalTo(fromBase(b2));
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
	public boolean equalTo(NucleotideIUPAC o) {
		switch (o) {
		case NaC:
			return false;
		case T:
		case U:
			return this == T || this == U;
		default:
			return this == o;
		}
	}

	/**
	 * Complementary code.
	 * 
	 * <p/>
	 * 
	 * The code that correspond to all nucleotides that are not included in this
	 * code. For example for {@link #A} would be {@link #B} (i.e. {@link #C},
	 * {@link #G} or {@link #T}).
	 * 
	 * <p/>
	 * There are a few special cases worth considering:
	 * 
	 * <ol>
	 * <li> <code> {@link #N}.complement() == {@link #NaC} </code></li>
	 * <li> <code> {@link #GAP}.complement() == {@link #NaC} </code></li>
	 * <li> <code> {@link #NaC}.complement() == {@link #NaC} </code></li>
	 * </ol>
	 * 
	 * Due to (1) and (2) complement is not invertible on {@link #N} nor on
	 * {@link #GAP}.
	 * 
	 * @param isDNA
	 *            indicates that it should return {@link #T} over {@link #U} as
	 *            a complement of {@code #V} ({@code true}) or <i>vice-versa</i>
	 *            ({@code false})
	 * 
	 * @return never {@code null}, but possible {@code #NaC} to indicate
	 *         non-sense complements.
	 */
	public NucleotideIUPAC complement(boolean isDNA) {
		switch (this) {
		case NaC:
		case N:
		case GAP:
			return NaC;
		case V:
			return isDNA ? T : U;
		default:
			return MASK_TO_CODE[(~this.maskValue) & N.maskValue];
		}
	}

	/**
	 * Complementary code assuming you are working on DNA.
	 * 
	 * Equivalent to
	 * <code> {@link #complement(boolean) complement(true)} </code>
	 * 
	 * @see #complement(boolean)
	 */
	public NucleotideIUPAC complement() {
		return this.complement(true);

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
