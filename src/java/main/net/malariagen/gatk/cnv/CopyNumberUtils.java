package net.malariagen.gatk.cnv;

import org.broadinstitute.sting.utils.GenomeLoc;

public class CopyNumberUtils {


	
	public static CopyNumberStatistics emptyStatistics() {
		return new CopyNumberStatistics() {
			
			private transient GenomeLoc location;
			@Override
			public GenomeLoc getLocation() {
				if (location != null)
					return location;
				return location = GenomeLoc.UNMAPPED;
			}
			
			@Override
			public int getReadStarts(int offset) {
				throw new IllegalArgumentException("offset out of range");
			}

			@Override
			public int[] getReadStarts() {
				return new int[0];
			}
		};
	}
}
