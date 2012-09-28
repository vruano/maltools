package net.malariagen.gatk.coverage;

import static net.malariagen.gatk.coverage.LocusCategory.*;
import net.malariagen.gatk.coverage.CoverageBiasWalker.GroupBy;
import net.malariagen.gatk.math.IntegerCountersIncrement;
import net.malariagen.utils.Nucleotide;

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
public class MaDPPcDistributionWalker extends IntegerStatDistributionWalker {
	
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
		int[] allNucDepth = new int[4];
		int[][] perGroupNucDepth = new int[result.groupValues.length][4];
		if (groupBy != GroupBy.NONE){
			for (PileupElement pe : pileup) {
				int baseIndex = Nucleotide.fromByte(pe.getBase()).ordinal();
				if (baseIndex <= Nucleotide.T.ordinal()) allNucDepth[baseIndex]++;
				if (groupBy.implies(GroupBy.RG)) {
					String readGroupName = pe.getRead().getReadGroup().getId();
					int readGroupIndex = groupIndices.get(readGroupName);
					if (baseIndex <= Nucleotide.T.ordinal()) perGroupNucDepth[readGroupIndex][baseIndex]++;
				}
				if (groupBy.implies(GroupBy.SM)) {
					String sampleName = pe.getRead().getReadGroup().getSample();
					int sampleIndex = groupIndices.get(sampleName);
					if (baseIndex <= Nucleotide.T.ordinal()) perGroupNucDepth[sampleIndex][baseIndex]++;
				}
			}
			for (int i = 0; i < perGroupNucDepth.length; i++) {
				int mjIndex = 0;
				int total = perGroupNucDepth[i][0];
				for (int j = 1; j < 4; j++) {
					int jDepth = perGroupNucDepth[i][j];
					total += jDepth;
					if (jDepth > perGroupNucDepth[i][mjIndex])
						mjIndex = j;
				}
				result.groupValues[i] = total == 0 ? 0 : (int) Math.round(100.0 * perGroupNucDepth[i][mjIndex]/(double)total);
			}
		}
		int mjIndex = 0;
		for (int j = 1; j < 4; j++) 
			if (allNucDepth[j] > allNucDepth[mjIndex])
				mjIndex = j;
		result.depth = result.depth == 0 ? 0 : (int) Math.round(100.0 * allNucDepth[mjIndex]/ (double) result.depth);
		return result;
	}

}
