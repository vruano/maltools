package net.malariagen.gatk.walker;

import net.sf.samtools.SAMRecord;

public enum FragmentFlag {

	FIRST_FACE_AWAY, LAST_FACE_AWAY, FIRST_UNMAPPED, LAST_UNMAPPED,
	FIRST_LOW_MQ, LAST_LOW_MQ, LARGE_FRAGMENT, CROSS_CHROMOSOME;
	
	public static int flagsOf(FragmentRecord fragment, int lowMqThr, int maxLength) {
		
		int result = 0;
		SAMRecord first = fragment.getFirstRead();
		SAMRecord last = fragment.getLastRead();
		boolean bothMapped = true;
		if (first.getReadUnmappedFlag()) { result |= 1 << FIRST_UNMAPPED.ordinal(); bothMapped = false; }
		if (last.getReadUnmappedFlag()) { result |= 1 << LAST_UNMAPPED.ordinal(); bothMapped = false; }
		if (lowMqThr > 0) {
			if (!first.getReadUnmappedFlag() && first.getMappingQuality() < lowMqThr) result |= 1 << FIRST_LOW_MQ.ordinal();
			if (!last.getReadUnmappedFlag() && last.getMappingQuality() < lowMqThr) result |= 1 << LAST_LOW_MQ.ordinal();
		}
		if (bothMapped) {
			SAMRecord refFirst = fragment.getReferenceFirstRead();
			SAMRecord refLast = fragment.getReferenceLastRead();
			if (refFirst.getReadNegativeStrandFlag()) 
				result |= 1 << (refFirst == first ? FIRST_FACE_AWAY.ordinal() : LAST_FACE_AWAY.ordinal());
			if (!refLast.getReadNegativeStrandFlag())
				result |= 1 << (refLast == last ? LAST_FACE_AWAY.ordinal() : FIRST_FACE_AWAY.ordinal());
			if (!refFirst.getReferenceName().equals(refLast.getReferenceName()))
				result |= 1 << CROSS_CHROMOSOME.ordinal();
			else if (fragment.getLength() > maxLength)
				result |= 1 << LARGE_FRAGMENT.ordinal();
		}
		return result;
	}
}
