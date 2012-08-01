package net.malariagen.gatk.filters;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.filters.ReadFilter;

public class UnmappedMateFilter extends ReadFilter {

	@Override
	public boolean filterOut(SAMRecord arg0) {
		if (!arg0.getReadPairedFlag()) return false;
		if (arg0.getMateUnmappedFlag()) return true;
		if (arg0.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) return true;
		if (arg0.getMateReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) return true;
		if (arg0.getMateAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) return true;
		return false;
	}

}
