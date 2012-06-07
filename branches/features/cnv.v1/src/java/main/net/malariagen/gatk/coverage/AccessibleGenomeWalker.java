package net.malariagen.gatk.coverage;

import net.malariagen.gatk.math.IntegerCounterSet;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

@ReadFilters(BadMateFilter.class)
@PartitionBy(PartitionType.LOCUS)
@By(DataSource.REFERENCE)
public class AccessibleGenomeWalker extends LocusWalker<LocusAccessibility, AccessibilityStatistics> implements
TreeReducible<IntegerCounterSet>{

	@Override
	public IntegerCounterSet treeReduce(IntegerCounterSet lhs,
			IntegerCounterSet rhs) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public LocusAccessibility map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public AccessibilityStatistics reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public AccessibilityStatistics reduce(LocusAccessibility value,
			AccessibilityStatistics sum) {
		// TODO Auto-generated method stub
		return null;
	}

}
