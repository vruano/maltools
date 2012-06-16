package net.malariagen.gatk.coverage;

import java.util.Collection;

import net.sf.samtools.SAMReadGroupRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;


public class FragmentLengthWalker extends ReadWalker<GATKSAMRecord,FragmentLengths> implements TreeReducible<FragmentLengths> {


	@Argument(shortName="maxfl", fullName="max_fragment_length", doc = "maximum admissible fragment length; mapped read pairs on the same cromosome that have a larger fragment size would be discarded", required = false)
	public int maxLength = 10000;
	
	@Override
	public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read,
			ReadMetaDataTracker metaDataTracker) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public FragmentLengths reduceInit() {
		Collection<? extends Sample> samples = getToolkit().getSampleDB().getSamples();
		SAMDataSource rds = getToolkit().getReadsDataSource();
		Collection<? extends SAMReadGroupRecord> rgs = 
		rds.getHeader().getReadGroups();
		return FragmentLengths.create(samples, rgs, maxLength);
	}

	@Override
	public FragmentLengths reduce(
			GATKSAMRecord value,
			FragmentLengths sum) {
		return FragmentLengths.add(value,sum);
	}

	@Override
	public FragmentLengths treeReduce(FragmentLengths lhs, FragmentLengths rhs) {
		return FragmentLengths.merge(lhs,rhs);
	}

}
