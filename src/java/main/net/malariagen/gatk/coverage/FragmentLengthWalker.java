package net.malariagen.gatk.coverage;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;


public class FragmentLengthWalker extends ReadWalker<GATKSAMRecord,FragmentLengths> {


	@Override
	public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read,
			ReadMetaDataTracker metaDataTracker) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public FragmentLengths reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public FragmentLengths reduce(
			GATKSAMRecord value,
			FragmentLengths sum) {
		// TODO Auto-generated method stub
		return null;
	}

}
