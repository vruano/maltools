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

import org.broadinstitute.sting.gatk.walkers.PartitionType;

@ReadFilters({ UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class })
@PartitionBy(PartitionType.LOCUS)
@By(DataSource.REFERENCE)
@Downsample(by = DownsampleType.NONE)
@Requires({ DataSource.READS, DataSource.REFERENCE_ORDERED_DATA })
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
public abstract class IntegerStatDistributionWalker extends
		LocusWalker<IntegerCountersIncrement, IntegerCounterSet> implements
		TreeReducible<IntegerCounterSet> {
	
	org.broadinstitute.sting.gatk.DownsampleType mx;
	@Output(doc = "File where to store the coverage distribution in JSON format", shortName = "o")
	@Gather(IntegerDistributionSetGatherer.class)
	private PrintWriter output;

	@Argument(fullName = "groupBy", shortName = "groupBy", doc = "Whether to collect statistics for each sample separatelly as well")
	public GroupBy groupBy = GroupBy.NONE;

	@Argument(fullName = "minimumBaseQuality", shortName = "mbq", doc = "Minimum quality for a base to be considered in coverage counting, -1 to indicate no minimum (default)", required = false)
	public int minimumBaseQuality = -1;

	@Argument(fullName = "minimumMappingQuality", shortName = "mmq", doc = "Minimum mapping quality for a read (thus its based) to be considered in coverage counting, -1 to indicate no minimum (default)", required = false)
	public int minimumMappingQuality = -1;

	@Argument(fullName = "minimumBAQ", shortName = "mbaq", doc = "Minimum BAQ for a base to be considered in coverage counting, -1 to indicate no minimum (default)", required = false)
	public int minimumBAQ = -1;

	@Argument(fullName = "excludeAmbigousRef", shortName = "exclAmbRef", doc = "Indicates whether totally degenerated ambiguous reference sites should be excluded", required = true)
	public boolean excludeAmbigousRef = false;

	@Argument(fullName = "excludeAmbigousCall", shortName = "exclAmbCall", doc = "Indicates whether totally degenerated ambiguous read-bases should be excluded", required = false)
	public boolean excludeAmbigousBase = false;

	@Argument(fullName = "includeReadDeletions", shortName = "inclRdDel", doc = "Consider also deletions in reads as covering the site", required = false)
	public boolean includeReadDeletions = false;

	@Argument(fullName = "features", shortName = "features", doc = "ROD data binding indicating site category features", required = false)
	protected RodBinding<GFFFeature> features = null;

	private String[] groupNames;
	private int groupCount;
	protected Map<String, Integer> groupIndices;

	private String[] sequenceNames;
	protected Map<String, Integer> sequenceIndices;

	protected PileupElementFilter pileupFilter;

	protected Set<IntegerCounterSet> counterSets;
	protected ThreadLocal<IntegerCounterSet> counterSetByThread;

	protected ConcurrentPool<IntegerCountersIncrement> incPool;

	private BAQ baqHMM;
	private IndexedFastaSequenceFile reference;
	private CalculationMode cmode;
	private QualityMode qmode;

	private ReadGroupDB readGroupDb;

	@Override
	public abstract IntegerCountersIncrement map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context);
	
	@Override
	public IntegerCounterSet reduceInit() {
		IntegerCounterSet result = counterSetByThread.get();
		if (result == null) {
			synchronized (counterSets) {
				counterSetByThread.set(result = new IntegerCounterSet(
						groupNames, sequenceNames));
				counterSets.add(result);
			}
		}
		return result;
	}

	@Override
	public void initialize() {
		groupInitialize();
		sequenceNamesInitialize();
		initializeCoverageCounterSets();
		initializeBAQCalculatingEngine();
		pileupFilter = initializePileupFilter(includeReadDeletions,
				excludeAmbigousBase, minimumBaseQuality, minimumMappingQuality, minimumBAQ,
				baqHMM, reference, cmode, qmode);
	}

	private void groupInitialize() {
		readGroupDb = new ReadGroupDB(this.getToolkit());
		Set<String> nameSet = new LinkedHashSet<String>();

		if (groupBy.implies(GroupBy.SM))
			nameSet.addAll(readGroupDb.getSampleIDs());
		if (groupBy.implies(GroupBy.RG))
			nameSet.addAll(readGroupDb.getReadGroupIDs());

		groupCount = nameSet.size();
		groupNames = nameSet.toArray(new String[groupCount]);
		if (groupCount == 0)
			groupIndices = Collections.emptyMap();
		else {
			int nextIndex = 0;
			groupIndices = new LinkedHashMap<String, Integer>();
			for (String groupName : groupNames)
				groupIndices.put(groupName, nextIndex++);
		}

	}

	private void initializeBAQCalculatingEngine() {
		baqHMM = new BAQ();
		reference = getToolkit().getReferenceDataSource().getReference();
		qmode = getToolkit().getWalkerBAQQualityMode();
		cmode = getToolkit().getArguments().BAQMode;
	}

	static PileupElementFilter initializePileupFilter(
			final boolean includeReadDeletions,
			final boolean excludeAmbigousBase, final int minimumBaseQuality,
			final int minimumMappingQuality, final int minimumBAQ,
			final BAQ baqHMM, final IndexedFastaSequenceFile reference,
			final CalculationMode cmode, final QualityMode qmode) {
		return new PileupElementFilter() {
			public boolean allow(PileupElement pe) {
				if (!includeReadDeletions
						&& pe.getBaseIndex() == BaseUtils.DELETION_INDEX)
					return false;
				if (excludeAmbigousBase
						&& !BaseUtils.isRegularBase(pe.getBase()))
					return false;
				if (minimumBaseQuality >= 0
						&& pe.getQual() < minimumBaseQuality)
					return false;
				if (minimumMappingQuality >= 0
						&& pe.getMappingQual() < minimumMappingQuality)
					return false;
				if (minimumBAQ >= 0
						&& baqFor(pe.getRead(), pe.getOffset(), baqHMM,
								reference, cmode, qmode) < minimumBAQ) {
					return false;
				}
				return true;
			}
		};
	}

	private static int baqFor(SAMRecord read, int offset, BAQ baqHMM,
			IndexedFastaSequenceFile reference, CalculationMode cmode,
			QualityMode qmode) {
		baqHMM.baqRead(
				read,
				reference,
				cmode == CalculationMode.RECALCULATE ? (BAQ.hasBAQTag(read) ? CalculationMode.CALCULATE_AS_NECESSARY
						: cmode)
						: cmode, qmode);
		return BAQ.calcBAQFromTag(read, offset, false);
	}

	private void initializeCoverageCounterSets() {
		counterSets = new HashSet<IntegerCounterSet>(20);
		counterSetByThread = new ThreadLocal<IntegerCounterSet>();
		incPool = new ConcurrentPool<IntegerCountersIncrement>(
				new Factory<IntegerCountersIncrement>() {
					public IntegerCountersIncrement newInstance() {
						IntegerCountersIncrement result = new IntegerCountersIncrement();
						result.groupValues = groupCount > 0 ? new int[groupCount]
								: null;
						return result;
					}
				}) {

			@Override
			public IntegerCountersIncrement borrow() {
				IntegerCountersIncrement result = super.borrow();
				result.clear();
				return result;
			}

		};

	}

	private void sequenceNamesInitialize() {
		SAMSequenceDictionary sd = getToolkit().getReferenceDataSource()
				.getReference().getSequenceDictionary();
		List<SAMSequenceRecord> sequences = sd.getSequences();
		sequenceNames = new String[sequences.size()];
		sequenceIndices = new HashMap<String, Integer>(sequences.size());
		int nextIndex = 0;
		for (SAMSequenceRecord s : sequences) {
			String name = s.getSequenceName();
			sequenceIndices.put(name, nextIndex);
			sequenceNames[nextIndex++] = name;
		}
	}

	@Override
	public IntegerCounterSet reduce(IntegerCountersIncrement cci,
			IntegerCounterSet ccs) {
		if (cci == null)
			return ccs;
		if (groupBy != GroupBy.NONE)
			ccs.addSampleValues(cci.groupValues, cci.categories, cci.sequence);
		else
			ccs.addAllValue(cci.depth, cci.categories, cci.sequence);
		incPool.restore(cci);
		return ccs;
	}

	@Override
	public IntegerCounterSet treeReduce(IntegerCounterSet lhs,
			IntegerCounterSet rhs) {
		// The per-thread reuse takes care of merging content
		return lhs;
	}

	@Override
	public void onTraversalDone(IntegerCounterSet sum) {
		IntegerCounterSet realSum = new IntegerCounterSet(groupNames,
				sequenceNames);
		for (IntegerCounterSet threadCCS : counterSets)
			realSum.applyCounterSet(threadCCS);
		logReport(realSum);

		try {
			realSum.toCoverageDistributionSet().write(output);
		} catch (IOException e) {
			throw new CoverageDistributionException(
					"could not write coverage distribution file", e);
		}
	}

	private void logReport(IntegerCounterSet sum) {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		sum.printReport(pw);
		pw.close();
		String report = sw.toString();
		String[] reportLines = report.split("\\n");
		for (String line : reportLines)
			logger.info(line);
	}

}