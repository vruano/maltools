package net.malariagen.gatk.uniqueness;

import org.broad.tribble.Feature;

public class UQNFeature implements Feature {

	private String sequence;
	
	private int position;
	
	private int score;
	
	public String getSequence() {
		return sequence;
	}
	
	public void setSequence(String seq) {
		sequence = seq;
	}
	
	public void setPosition(int pos) {
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

	public int getScore() {
		return score;
	}
	
	public void setScore(int score) {
		this.score = score;
	}
}
