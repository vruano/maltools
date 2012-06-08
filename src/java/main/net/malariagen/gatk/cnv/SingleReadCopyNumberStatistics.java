package net.malariagen.gatk.cnv;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.utils.GenomeLoc;

public class SingleReadCopyNumberStatistics implements CopyNumberStatistics {

	private GenomeLoc loc;
	
	public SingleReadCopyNumberStatistics(GenomeLoc loc) {
		this.loc = loc;
		if (loc.size() != 1)
			throw new IllegalArgumentException("the location provided must have size 1");
	}

	@Override
	public GenomeLoc getLocation() {
		return loc;
	}

	@Override
	public int getReadStarts(int offset) {
		if (offset == 0)
			return 1;
		throw new IllegalArgumentException("offset out of range");
	}

	@Override
	public int[] getReadStarts() {
		return new int[] { 1 };
	}

}
