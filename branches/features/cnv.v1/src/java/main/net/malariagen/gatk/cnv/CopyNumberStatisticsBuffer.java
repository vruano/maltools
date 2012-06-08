package net.malariagen.gatk.cnv;

import java.util.HashMap;
import java.util.Map;

import org.broadinstitute.sting.utils.GenomeLoc;

import scala.actors.threadpool.Arrays;

public class CopyNumberStatisticsBuffer {
	

	public Map<String,ContigBuffer> contigBuffers = new HashMap<String,ContigBuffer>(20);
	
	public CopyNumberStatisticsBuffer add(CopyNumberStatistics css) {
		
		GenomeLoc loc = css.getLocation();
		if (GenomeLoc.isUnmapped(loc)) return this;

		ContigBuffer contigBuffer = getContigBuffer(loc.getContig());
		contigBuffer.add(css);
		return this;
	}

	protected CopyNumberStatistics get(GenomeLoc loc) {
		String contig = loc.getContig();
		ContigBuffer counter = contigBuffers.get(contig);
		if (loc == null)
			return CopyNumberUtils.emptyStatistics();
		return counter.get(loc);
	}
	
	
	private ContigBuffer getContigBuffer(String contig) {
		ContigBuffer result = contigBuffers.get(contig);
		if (result == null)
			contigBuffers.put(contig,result = new ContigBuffer());
		return result;
	}

	protected class ContigBuffer {

		private int readStarts[];
		
		public void add(CopyNumberStatistics css) {
			
			GenomeLoc loc = css.getLocation();
			int start = loc.getStart();
			int stop = loc.getStop();
			ensureCapacity(start,stop);
			int[] cssRs = css.getReadStarts();
			for (int i = 0; i < cssRs.length; i++)
				readStarts[start - 1 + i] += cssRs[i];
			
		}

		public CopyNumberStatistics get(GenomeLoc loc) {
			int start = loc.getStart();
			int stop = loc.getStop();
			int[] rs = Arrays.copyOfRange(readStarts, start - 1, stop);
			return new SimpleCopyNumberStatistics(loc,rs);
		}


		private void ensureCapacity(int start, int stop) {
			if (readStarts == null)
				readStarts = new int[stop];
			else if (readStarts.length < stop)
				readStarts = Arrays.copyOf(readStarts, stop);
		}
	}
	
}
