package net.malariagen.gatk.math;

import java.util.Arrays;

import net.malariagen.gatk.coverage.LocusCategory;



public class IntegerCountersIncrement {

	public IntegerCountersIncrement () {}
	
	//total depth at that site.
	public int depth;

	//applicable counter categories.
	public byte categories;

	//sequence name.
	public int sequence;

	//depth per sample. could be null if we are not recording per sample counters.
	public int[] groupValues;
	
	//generates the mask given a set of locus categories
	public static byte categoryMask(LocusCategory ... categories) {
		byte result = 0;
		for (LocusCategory c : categories) {
			int ord = c.ordinal();
			result |= 1 << ord;
		}
		return result;
	}

	public void addCategory(LocusCategory c) {
		categories |= 1 << c.ordinal();
	}

	public void clear() {
		if (groupValues != null)
			Arrays.fill(groupValues, 0);
		depth = 0;
		categories = 0;
		sequence = -1;
	}
	
}
