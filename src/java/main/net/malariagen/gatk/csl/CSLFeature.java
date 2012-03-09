package net.malariagen.gatk.csl;

import org.broad.tribble.Feature;

public class CSLFeature implements Feature {

	private String sequence;
	
	private int position;
	
	private byte reference;
	
	public byte getReference() {
		return reference;
	}

	void setReference(byte reference) {
		this.reference = reference;
	}

	public byte getAlternative() {
		return alternative;
	}

	 void setAlternative(byte alternative) {
		this.alternative = alternative;
	}

	private byte alternative;
	
	
	public String getSequence() {
		return sequence;
	}
	
	 void setSequence(String seq) {
		sequence = seq;
	}
	
	void setPosition(int pos) {
		position = pos; 
	}
	
	@Override
	public String getChr() {
		return sequence;
	}

	@Override
	public int getEnd() {
		return position;
	}

	@Override
	public int getStart() {
		return position;
	}


}
