package net.malariagen.gatk.cnv;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

public class CopyNumberWalker extends ReadWalker<CopyNumberStatistics,CopyNumberStatisticsBuffer> {


	@Override
	public CopyNumberStatisticsBuffer reduceInit() {
		return new CopyNumberStatisticsBuffer();
	}

	@Override
	public CopyNumberStatisticsBuffer reduce(CopyNumberStatistics value,
			CopyNumberStatisticsBuffer sum) {
		
		return sum.add(value);
	}

	@Override
	public CopyNumberStatistics map(ReferenceContext ref, SAMRecord read,
			ReadMetaDataTracker metaDataTracker) {
		GenomeLoc loc = ref.getLocus();
		loc = ref.getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart(), loc.getStart());
		return new SingleReadCopyNumberStatistics(loc);
	}
	
	@Override
	public void onTraversalDone(CopyNumberStatisticsBuffer buffer) {
		for (GenomeLoc contigLocation : buffer.contigLocationSet()) {
			analyzeContig(buffer.get(contigLocation));
		}
	}

	private void analyzeContig(CopyNumberStatistics copyNumberStatistics) {
		// TODO Auto-generated method stub
		
	}

}
