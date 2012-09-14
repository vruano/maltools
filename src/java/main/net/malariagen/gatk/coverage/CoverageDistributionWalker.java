package net.malariagen.gatk.coverage;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static net.malariagen.gatk.coverage.LocusCategory.*;
import net.malariagen.gatk.coverage.CoverageBiasWalker.GroupBy;
import net.malariagen.gatk.gff.GFFFeature;
import net.malariagen.gatk.math.IntegerCounterSet;
import net.malariagen.gatk.math.IntegerCountersIncrement;
import net.malariagen.gatk.math.IntegerDistributionSetGatherer;
import net.malariagen.gatk.utils.ReadGroupDB;
import net.malariagen.utils.ConcurrentPool;
import net.malariagen.utils.Factory;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Gather;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
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
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode;
import org.broadinstitute.sting.utils.baq.BAQ.QualityMode;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElementFilter;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import org.broadinstitute.sting.gatk.walkers.PartitionType;

@ReadFilters({ UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class })
@PartitionBy(PartitionType.LOCUS)
@By(DataSource.REFERENCE)
@Downsample(by = DownsampleType.NONE)
@Requires({ DataSource.READS, DataSource.REFERENCE_ORDERED_DATA })
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
public class CoverageDistributionWalker extends IntegerStatDistributionWalker {
	

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
		if (groupBy != GroupBy.NONE) 
			for (PileupElement pe : pileup) {
				if (groupBy.implies(GroupBy.RG)) {
					String readGroupName = pe.getRead().getReadGroup().getId();
					int readGroupIndex = groupIndices.get(readGroupName);
					result.groupValues[readGroupIndex]++;
				}
				if (groupBy.implies(GroupBy.SM)) {
					String sampleName = pe.getRead().getReadGroup().getSample();
					int sampleIndex = groupIndices.get(sampleName);
					result.groupValues[sampleIndex]++;
				}
			}
		return result;
	}

}
