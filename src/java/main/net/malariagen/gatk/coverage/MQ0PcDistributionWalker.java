package net.malariagen.gatk.coverage;

import static net.malariagen.gatk.coverage.LocusCategory.*;
import net.malariagen.gatk.coverage.CoverageBiasWalker.GroupBy;
import net.malariagen.gatk.math.IntegerCountersIncrement;

import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Requires;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import org.broadinstitute.sting.gatk.walkers.PartitionType;

@ReadFilters({ UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class })
@PartitionBy(PartitionType.LOCUS)
@By(DataSource.REFERENCE)
@Downsample(by = DownsampleType.NONE)
@Requires({ DataSource.READS, DataSource.REFERENCE_ORDERED_DATA })
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
public class MQ0PcDistributionWalker extends IntegerStatDistributionWalker {
	

	@Override
	public IntegerCountersIncrement map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context) {

		if (excludeAmbigousRef && !BaseUtils.isRegularBase(ref.getBase()))
			return null;
		IntegerCountersIncrement result = incPool.borrow();
		result.categories = categoryMask(tracker, features);
		result.sequence = sequenceIndices.get(ref.getLocus().getContig());
		ReadBackedPileup pileup = context.getBasePileup().getFilteredPileup(
				pileupFilter);
		result.depth = context.size();
		int totalMQ0 = 0;
		if (groupBy != GroupBy.NONE){
			int[] mq0 = new int[result.groupValues.length];
			for (PileupElement pe : pileup) {
				if (groupBy.implies(GroupBy.RG)) {
					String readGroupName = pe.getRead().getReadGroup().getId();
					int readGroupIndex = groupIndices.get(readGroupName);
					result.groupValues[readGroupIndex]++;
					if (pe.getMappingQual() == 0)
						mq0[readGroupIndex]++;
				}
				if (groupBy.implies(GroupBy.SM)) {
					String sampleName = pe.getRead().getReadGroup().getSample();
					int sampleIndex = groupIndices.get(sampleName);
					result.groupValues[sampleIndex]++;
					if (pe.getMappingQual() == 0)
						mq0[sampleIndex]++;
				}
			}
			for (int i = 0; i < mq0.length; i++)
				result.groupValues[i] = (result.groupValues[i]== 0) ? 0 :  (mq0[i] * 100) / result.groupValues[i];
		}
		result.depth = result.depth == 0 ? 0 : (totalMQ0 * 100) / result.depth;
		return result;
	}

}
