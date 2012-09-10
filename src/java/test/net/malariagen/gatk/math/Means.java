package net.malariagen.gatk.math;

import org.apache.commons.math.stat.descriptive.rank.Median;

public final class Means {
	
	public static double rms(int ... values) {
		long sqSum = 0;
		for (int i : values) 
			sqSum += i * i;
		return Math.sqrt(((double) sqSum)/values.length);
	}
	
	public static double rms(int[] values, int fromIndex, int toIndex) {
		if (fromIndex < 0 || fromIndex >= values.length)
			throw new IllegalArgumentException("fromIndex must be from 0 to the input length - 1");
		if (toIndex <0 || toIndex >= values.length)
			throw new IllegalArgumentException("toIndex must be from 0 to the input length");
		if (fromIndex == 0 && toIndex == values.length)
			return rms(values);
		else {
			long sqSum = 0;
			for (int i = fromIndex; i < toIndex; i++)
				sqSum += values[i] * values[i];
			return Math.sqrt(((double) sqSum)/values.length);
		}	
	}
	

}
