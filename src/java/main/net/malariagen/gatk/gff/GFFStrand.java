package net.malariagen.gatk.gff;

/**
 * Gff feature strand.
 * @author vrr
 */
public enum GFFStrand {
	PLUS('+'), MINUS('-'), NA('.');

	private char charCode;

	
	// Creates the strand given its single char code.
	GFFStrand(char code) {
		this.charCode = code;
	}

	/**
	 * Parses the strand from a character sequence
	 * 
	 * @param cs
	 *            the character sequence to parse into an strand.
	 * @throws NullPointerException
	 *             if {@code cs} is {@code null}.
	 * @throws IllegalArgumentException
	 *             if the input char sequence does not have length 1.
	 * @return never {@code null} but {@link #NA} for non-applicable/empty ('.')
	 *         strands.
	 */
	public static GFFStrand fromCharSequence(CharSequence cs) {
		if (cs.length() != 1)
			throw new IllegalArgumentException(
					"input char sequence must have length 1");
		return fromChar(cs.charAt(0));
	}

	/**
	 * Transforms a char code into the corresponding strand.
	 * 
	 * @param code
	 * @throws IllegalArgumentException
	 *             if {@code code} is not a valid strand code.
	 * @return never {@code null}
	 */
	public static GFFStrand fromChar(char code) {
		switch (code) {
		case '+':
			return PLUS;
		case '-':
			return MINUS;
		case '.':
			return NA;
		default:
			throw new IllegalArgumentException("Illegal strand code");
		}
	}

	/**
	 * Returns the single character code for this strand.
	 */
	public char toChar() {
		return charCode;
	}

	@Override
	public String toString() {
		return "" + charCode;
	}
}
