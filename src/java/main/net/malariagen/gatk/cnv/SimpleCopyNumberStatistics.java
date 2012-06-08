package net.malariagen.gatk.cnv;

import org.broadinstitute.sting.utils.GenomeLoc;

public class SimpleCopyNumberStatistics implements CopyNumberStatistics {

	private GenomeLoc loc;
	private int[] rs;

	public SimpleCopyNumberStatistics(GenomeLoc loc, int[] rs) {
		this.loc = loc;
		this.rs = rs;
	}

	@Override
	public GenomeLoc getLocation() {
		return loc;
	}

	@Override
	public int[] getReadStarts() {
		return rs;
	}

	@Override
	public int getReadStarts(int offset) {
		return rs[offset];
	}

}
