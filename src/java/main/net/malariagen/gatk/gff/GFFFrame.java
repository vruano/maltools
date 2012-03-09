package net.malariagen.gatk.gff;

/**
 * Gff feature frame offset.
 * @author vrr
 */
public enum GFFFrame {
	ZERO('0'), ONE('1'), TWO('2'), NA('.');
	
	private char charCode;
	
	GFFFrame(char code) {
		charCode = code;
	}
	
	public int toOffset() {
		return this.ordinal();
	}
		
	

	/**
	 * Parses the frame from a character sequence
	 * 
	 * @param cs
	 *            the character sequence to parse into a frame.
	 * @throws NullPointerException
	 *             if {@code cs} is {@code null}.
	 * @throws IllegalArgumentException
	 *             if the input char sequence does not have length 1.
	 * @return never {@code null} but {@link #NA} for non-applicable/empty ('.')
	 *         strands.
	 */
	public static GFFFrame fromCharSequence(CharSequence cs) {
		if (cs.length() != 1)
			throw new IllegalArgumentException(
					"input char sequence must have length 1");
		return fromChar(cs.charAt(0));
	}
	
	

	/**
	 * Transforms a char code into the corresponding frame.
	 * 
	 * @param code
	 * @throws IllegalArgumentException
	 *             if {@code code} is not a valid strand code.
	 * @return never {@code null}
	 */
	public static GFFFrame fromChar(char code) {
		switch (code) {
		case '0' : return ZERO;
		case '1' : return ONE;
		case '2' : return TWO;
		case '.' : return NA;
		default: 
			throw new IllegalArgumentException("char-code '" + code + "' not valid for a gff frame");
		}
	}
	
	public GFFFrame fromOffset(int offset) {
		switch (offset) {
		case 0 : return ZERO;
		case 1 : return ONE;
		case 2 : return TWO;
		default : 
			throw new IllegalArgumentException("offset '" + offset + "' not valid for a gff frame");
		}
	}
	
	public char toChar() {
		return charCode;
	}
	
	@Override
	public String toString() {
		return "" + charCode;
	}
	
}
