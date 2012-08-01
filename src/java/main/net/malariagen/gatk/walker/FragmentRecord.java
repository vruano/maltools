package net.malariagen.gatk.walker;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class FragmentRecord {
	
	public GATKSAMRecord getFirstRead() {
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
	
	GATKSAMRecord getFirstReadReported() {
		return swapped ? second : first;
	}

	private GATKSAMRecord first;
	private GATKSAMRecord second;
	private boolean swapped;
	private int fragmentLength;
	
	FragmentRecord(GATKSAMRecord r1, GATKSAMRecord r2, int fl) {
		if (r1.getFirstOfPairFlag()) {
		  first = r1;
		  second = r2;
		  swapped = false;
		}
		else {
		  first = r2;
		  second = r1;
		  swapped = true;
		}
		fragmentLength = fl;
	}
	
	public int getLength() {
		return fragmentLength;
	}
	
	public boolean getPositiveStrand() {
		return ! first.getReadNegativeStrandFlag();
	}
	
	
	/**
	 * Returns the read record of the fragment that appears first in the reference genome. 
	 * 
	 * If either read in unmapped, the result of this method is undetermined.
	 * 
	 * @return never {@code null}
	 */
	public GATKSAMRecord getReferenceFirstRead() {
		return first.getAlignmentStart() <= second.getAlignmentStart() ? first : second;
	}

	/**
	 * Returns the read record of the fragment that appears last in the reference genome. 
	 * 
	 * If either read in unmapped, the result of this method is undetermined.
	 * 
	 * @return never {@code null}
	 */
	public GATKSAMRecord getReferenceLastRead() {
		return first.getAlignmentStart() <= second.getAlignmentStart() ? second : first;
	}
	
	/**
	 * Returns the region expanding the fragment or at least the mapped reads that belong to it.
	 * <p/>
	 * If there is no mapped read it returns {@link GenomeLoc#UNMAPPED}
	 * <p/>
	 * 
	 * @param IllegalStateException 
	 *    if the pair map to different chromosomes, the result in undetermined.
	 * 
	 * @param parser the location parser to use in order to create the location.
	 * 
	 * @return never {@code null}. 
	 */
	public GenomeLoc getRegion(GenomeLocParser parser) {
		String contig = first.getReferenceName();
		if (first.getReadUnmappedFlag() && second.getReadUnmappedFlag())
			return GenomeLoc.UNMAPPED;
		else if (first.getReadUnmappedFlag())
			return parser.createGenomeLoc(second.getReferenceName(),second.getAlignmentStart(),second.getAlignmentStart());
		else if (second.getReadUnmappedFlag())
			return parser.createGenomeLoc(first.getReferenceName(),first.getAlignmentStart(),first.getAlignmentStart());
		int start = Math.min(first.getAlignmentStart(),second.getAlignmentStart());
		int stop = Math.max(first.getAlignmentEnd(),second.getAlignmentEnd());
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

	public boolean isTranslocation(long maxLength) {
		if (!isMapped()) return false;
		if (!this.first.getReferenceName().equals(this.second.getReferenceName()))
			return true;
		if (this.getLength() > maxLength)
			return true;
		return false;
	}
	
	public boolean isTranslocation() {
		return this.isTranslocation(Long.MAX_VALUE);
	}

	public boolean isMapped() {
		return ! first.getReadUnmappedFlag() && ! second.getReadUnmappedFlag();
	}
	
	public String toString() {
		return super.toString();
	}

	public GenomeLoc getFragmentStart(GenomeLocParser locParser) {
		if (locParser == null)
			throw new IllegalArgumentException("illegal locParser; must not be null");
		SAMRecord first = this.getFirstRead();
		if (first.getReadUnmappedFlag())
			return GenomeLoc.UNMAPPED;
		if (first.getReadNegativeStrandFlag())
			return locParser.createGenomeLoc(first.getReferenceName(),first.getAlignmentEnd(),first.getAlignmentEnd());
		else
			return locParser.createGenomeLoc(first.getReferenceName(),first.getAlignmentStart(),first.getAlignmentStart());
	}

	public GenomeLoc getReferenceStart(GenomeLocParser locParser) {
		if (locParser == null)
			throw new IllegalArgumentException("illegal locParser; must not be null");
		SAMRecord first = getReferenceFirstRead();
		if (first.getReadUnmappedFlag())
			return GenomeLoc.UNMAPPED;
		else 
			return locParser.createGenomeLoc(first.getReferenceName(),first.getAlignmentStart(),first.getAlignmentStart());
	}

	GATKSAMRecord getLastReadReported() {
		return swapped ? first : second;
	}
}
