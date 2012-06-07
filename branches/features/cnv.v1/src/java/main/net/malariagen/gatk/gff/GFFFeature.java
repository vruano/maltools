package net.malariagen.gatk.gff;

import org.broad.tribble.Feature;

public class GFFFeature implements Feature {

	public static final String UNKNOWN_SEQUENCE = "<unknown>";
	
	public static final String UNKNOWN_SOURCE = "<unknown>";
	
	public static final Double NO_SCORE = Double.NaN;

	public String getSequenceName() {
		return sequenceName;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}

	public GFFType getType() {
		return type;
	}

	public void setType(GFFType type) {
		this.type = type;
	}

	public String getTypeString() {
		return typeString;
	}

	/**
	 * Changes the feature type string. It also changes the type object accoordingly.
	 * 
	 * @param s the new feature type string
	 * @throws NullPointerException is {@code s} is {@code null}.
	 */
	public void setTypeString(String s) {
		if (s == null)
			throw new NullPointerException();
		if (s.equals(typeString)) 
			return;
		this.typeString = s;
		this.type = GFFType.fromCharSequence(s);
	}

	public Double getScore() {
		return score;
	}

	public void setScore(Double score) {
		this.score = score;
	}
	
	/**
	 * Parses in the score from a char-sequence.
	 * 
	 * @throws NumberFormatException if the score is not empty (.) and it is in an invalid number format.
	 */
	public void setScore(CharSequence cs) {
		if (cs.length() == 1 && cs.charAt(0) == '.') 
			setScore(NO_SCORE);
		else {
			setScore(Double.parseDouble(cs.toString()));
		}
	}
	
	public GFFStrand getStrand() {
		return strand;
	}

	public void setStrand(GFFStrand strand) {
		this.strand = strand;
	}

	public GFFFrame getFrame() {
		return frame;
	}

	public void setFrame(GFFFrame frame) {
		this.frame = frame;
	}

	public String getAttributeString() {
		return attributeString;
	}

	public void setAttributeString(String attributeString) {
		this.attributeString = attributeString;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	private String sequenceName = UNKNOWN_SEQUENCE;
	
	private String source = UNKNOWN_SOURCE;
	
	private GFFType type = GFFType.UNKNOWN;
	
	private String typeString = type.toString();
	
	private int start = 1;
	
	private int end = 1;
	
	private Double score = null; // no score
	
	private GFFStrand strand = GFFStrand.NA;
	
	private GFFFrame frame = GFFFrame.NA;
	
	private String attributeString = "";
	
	
	
	
	@Override
	public String getChr() {
		return sequenceName;
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getEnd() {
		return end;
	}

}
