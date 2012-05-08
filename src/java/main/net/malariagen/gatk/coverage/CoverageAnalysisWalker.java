package net.malariagen.gatk.coverage;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;



public class CoverageAnalysisWalker extends RodWalker<LocusCoverageAnalysis,IntervalCoverageAnalysis> {

	
	@Override
	public void initialize () {
		super.initialize();
	}
	
	@Override
	public LocusCoverageAnalysis map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public IntervalCoverageAnalysis reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public IntervalCoverageAnalysis reduce(LocusCoverageAnalysis value,
			IntervalCoverageAnalysis sum) {
		// TODO Auto-generated method stub
		return null;
	}

}
