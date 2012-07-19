package net.malariagen.gatk.coverage;

import java.util.List;

import net.sf.samtools.AlignmentBlock;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class FragmentRecord {
	
	public GATKSAMRecord getStartRead() {
		return first;
	}

	public void setFirst(GATKSAMRecord first) {
		this.first = first;
	}

	public GATKSAMRecord getLastRead() {
		return second;
	}

	public void setSecond(GATKSAMRecord second) {
		this.second = second;
	}

	private GATKSAMRecord second;
	private GATKSAMRecord first;
	private int fragmentLength;
	
	FragmentRecord(GATKSAMRecord r1, GATKSAMRecord r2, int fl) {
		if (r1.getFirstOfPairFlag()) {
		  first = r1;
		  second = r2;
		}
		else {
		  first = r2;
		  second = r1;
		}
		fragmentLength = fl;
	}
	
	public int getLength() {
		return fragmentLength;
	}
	
	public boolean getPositiveStrand() {
		return ! first.getReadNegativeStrandFlag();
	}
	
	public GATKSAMRecord getReferenceFirstRead() {
		return ! first.getReadNegativeStrandFlag() ? first : second;
	}

	public GATKSAMRecord getReferenceLastRead() {
		return first.getReadNegativeStrandFlag() ? first : second;
	}
	
	public GenomeLoc getRegion(GenomeLocParser parser) {
		String contig = first.getReferenceName();
		int start = getReferenceFirstRead().getAlignmentBlocks().get(0).getReferenceStart();
		GATKSAMRecord endRecord = getReferenceLastRead();
		int stop = lastReferencePosition(endRecord);
		GenomeLoc result = parser.createGenomeLoc(contig, start, stop); 
		return result;
	}

	private int lastReferencePosition(GATKSAMRecord endRecord) {
//		List<AlignmentBlock> abs = endRecord.getAlignmentBlocks();
//		AlignmentBlock ab = abs.get(abs.size() -1);
//		int stop = ab.getLength() + ab.getReferenceStart() - 1;
//		return stop;
		return endRecord.getAlignmentEnd();
	}
	
	public String getContig() {
		return first.getReferenceName();
	}
	
	public int getStart() {
		return getReferenceFirstRead().getAlignmentStart();
	}
	
	public int getStop() {
		return getReferenceLastRead().getAlignmentEnd();
	}
	
	public int getInsertStart() {
		return lastReferencePosition(getReferenceFirstRead()) + 1;
	}
	
	public int getInsertEnd() {
		return getReferenceLastRead().getAlignmentStart() - 1;			
	}
	
	public int getInsertLength() {
		return getInsertEnd() - getInsertStart() + 1;
	}
}
