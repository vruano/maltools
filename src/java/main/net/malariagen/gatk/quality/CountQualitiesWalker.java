package net.malariagen.gatk.quality;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static net.malariagen.gatk.coverage.LocusCategory.*;
import net.malariagen.gatk.math.IntegerCounterSet;
import net.malariagen.gatk.math.IntegerDistributionSetGatherer;
import net.malariagen.utils.ConcurrentPool;
import net.malariagen.utils.Factory;
import net.malariagen.utils.Pool;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Gather;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode;
import org.broadinstitute.sting.utils.baq.BAQ.QualityMode;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@Requires({ DataSource.READS, DataSource.REFERENCE_BASES })
@ReadFilters({ BadMateFilter.class, UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class })
@PartitionBy(PartitionType.INTERVAL)
@By(DataSource.READS)
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
public class CountQualitiesWalker extends
		ReadWalker<QualityCountersIncrement, QualityCounters> implements
		TreeReducible<QualityCounters> {

	@Output(doc = "File where to store the BAQ distribution in JSON format", fullName = "baqOutput", shortName = "baqOut", required = false)
	@Gather(IntegerDistributionSetGatherer.class)
	private PrintWriter baqOutput;

	@Output(doc = "File where to store the mapping quality distribution in JSON format", fullName = "mappingQualityOutput", shortName = "mqOut", required = false)
	@Gather(IntegerDistributionSetGatherer.class)
	private PrintWriter mappingQualityOutput;

	@Output(doc = "File where to store the base quality distribution in JSON format", fullName = "baseQualityOutput", shortName = "bqOut", required = false)
	@Gather(IntegerDistributionSetGatherer.class)
	private PrintWriter baseQualityOutput;

	@Argument(fullName = "minimumMappingQualityForBaseCounting", shortName = "minMq4Bc", doc = "minimum mapping quality for a read for its base qualities to be considered, none by default (-1)", required = false)
	private int minimumMappingQualityForBaseCounting = -1;

	@Argument(fullName = "excludeUnmappedBases", shortName = "exclUmB", doc = "whether to exclude bases that are not mapped to the reference (e.g. they are found within and read-insertion)", required = false)
	private boolean excludeUnmappedBases = false;

	@Argument(fullName = "excludeAmbigousBases", shortName = "exclAmB", doc = "Indicates whether totally degenerated ambiguous read-bases should be excluded", required = false)
	public boolean excludeAmbigousBases = false;

	@Argument(fullName = "countSamples", shortName = "cs", doc = "Whether to collect statistics for each sample separatelly as well", required = false)
	public boolean countSamples = true;

	private String[] sampleNames;
	private int sampleCount;
	private Map<String, Integer> sampleIndices;
	private String[] sequenceNames;

	protected Set<QualityCounters> allQualityCounters;
	protected ThreadLocal<QualityCounters> qualityCountersByThread;

	private Map<String, Integer> sequenceIndices;

	private Pool<QualityCountersIncrement> incPool = new ConcurrentPool<QualityCountersIncrement>(
			new Factory<QualityCountersIncrement>() {
				public QualityCountersIncrement newInstance() {
					QualityCountersIncrement result = new QualityCountersIncrement();
					return result;
				}

			});

	private BAQ baqHMM;
	private IndexedFastaSequenceFile reference;
	private CalculationMode cmode;
	private QualityMode qmode;

	

	@Override
	public QualityCountersIncrement map(ReferenceContext ref,
			GATKSAMRecord read, ReadMetaDataTracker tracker) {

		int mappingQuality = read.getMappingQuality();
		if (mappingQuality < minimumMappingQualityForBaseCounting
				&& mappingQualityOutput == null)
			return null;
		QualityCountersIncrement result = incPool.borrow();
		result.ensureCapacity(read.getReadLength());
		result.categories = categoryMask(tracker, read);
		result.sequence = sequenceIndices.get(ref.getLocus().getContig());
		result.sample = sampleIndices.get(read.getReadGroup().getSample());
		result.mappingQuality = mappingQuality;
		if (result.mappingQuality < minimumMappingQualityForBaseCounting
				|| (baseQualityOutput == null && baqOutput == null))
			result.baseQualityCount = result.baqCount = 0;
		else {
			byte[] baseQualities = read.getBaseQualities();
			if (baqOutput != null)
				baqHMM.baqRead(read, reference, cmode, qmode);
			byte[] baqs = BAQ.calcBAQFromTag(read, false, baqOutput == null);

			int nextIndex = 0;
			byte[] readBases = excludeAmbigousBases ? read.getReadBases() : null;
			if (this.excludeUnmappedBases)
				for (AlignmentBlock block : read.getAlignmentBlocks()) {
					int start = block.getReadStart() - 1;
					int end = start + block.getLength();
					for (int i = start; i < end; i++) {
						if (readBases != null
								&& !BaseUtils.isRegularBase(readBases[i]))
							continue;
						result.baseQualities[nextIndex] = baseQualities[i];
						result.baqs[nextIndex++] = baqs[i];
					}
				}
			else
				for (int i = 0; i < baseQualities.length; i++) {
					if (readBases != null
							&& !BaseUtils.isRegularBase(readBases[i]))
						continue;
					result.baseQualities[nextIndex] = baseQualities[i];
					result.baqs[nextIndex++] = baqs[i];
				}
			result.baseQualityCount = result.baqCount = nextIndex;
		}
		return result;
	}

	@Override
	public QualityCounters reduceInit() {
		QualityCounters result = qualityCountersByThread.get();
		if (result == null) {
			synchronized (allQualityCounters) {
				qualityCountersByThread.set(result = new QualityCounters(
						mappingQualityOutput != null,
						baseQualityOutput != null, baqOutput != null,
						sampleNames, sequenceNames));
				allQualityCounters.add(result);
			}
		}
		return result;
	}

	@Override
	public void initialize() {
		if (baqOutput == null && mappingQualityOutput == null
				&& baseQualityOutput == null)
			throw new RuntimeException(
					"at least one of quality output must be requested");
		sampleNamesInitialize();
		sequenceNamesInitialize();
		initializeBAQCalculatingEngine();
		initializeCoverageCounterSets();
	}

	private void initializeBAQCalculatingEngine() {
		baqHMM = new BAQ();
		reference = getToolkit().getReferenceDataSource().getReference();
		qmode = getToolkit().getWalkerBAQQualityMode();
		cmode = getToolkit().getArguments().BAQMode;
	}

	private void initializeCoverageCounterSets() {
		allQualityCounters = new HashSet<QualityCounters>(20);
		qualityCountersByThread = new ThreadLocal<QualityCounters>();
	}

	private void sequenceNamesInitialize() {
		SAMSequenceDictionary sd = getToolkit().getReferenceDataSource()
				.getReference().getSequenceDictionary();
		List<SAMSequenceRecord> sequences = sd.getSequences();
		sequenceNames = new String[sequences.size()];
		sequenceIndices = new HashMap<String, Integer>(sequenceNames.length);
		int nextIndex = 0;
		for (SAMSequenceRecord s : sequences) {
			String name = s.getSequenceName();
			sequenceIndices.put(name, nextIndex);
			sequenceNames[nextIndex++] = name;
		}
	}

	private void sampleNamesInitialize() {
		if (!countSamples) {
			sampleCount = 0;
			sampleNames = new String[0];
			sampleIndices = Collections.emptyMap();
		} else {
			LinkedList<String> allNames = new LinkedList<String>();
			for (Sample sample : getToolkit().getSampleDB().getSamples()) {
				allNames.add(sample.getID());
			}
			sampleCount = allNames.size();
			sampleNames = allNames.toArray(new String[sampleCount]);
			sampleIndices = new HashMap<String, Integer>(sampleCount);
			for (int i = 0; i < sampleCount; i++)
				sampleIndices.put(sampleNames[i], i);
		}
	}

	@Override
	public QualityCounters reduce(QualityCountersIncrement cci,
			QualityCounters ccs) {
		if (cci == null)
			return ccs;
		cci.applyTo(ccs, countSamples);
		incPool.restore(cci);
		return ccs;
	}

	@Override
	public QualityCounters treeReduce(QualityCounters lhs, QualityCounters rhs) {
		// lhs.applyCounterSet(rhs);
		return lhs;
	}

	@Override
	public void onTraversalDone(QualityCounters sum) {
		QualityCounters realSum;
		if (allQualityCounters.size() == 1)
			realSum = allQualityCounters.iterator().next();
		else {
			realSum = new QualityCounters(
			mappingQualityOutput != null, baseQualityOutput != null,
					baqOutput != null, sampleNames, sequenceNames);
			for (QualityCounters threadCCS : allQualityCounters)
				realSum.addFrom(threadCCS);
		}
		logReport(realSum);

		try {
			realSum.write(mappingQualityOutput, baseQualityOutput, baqOutput);
		} catch (IOException e) {
			throw new CountQualityException(
					"could not write coverage distribution file", e);
		}
	}

	private void logReport(QualityCounters qc) {
		if (qc.mappingQuality != null) {
			logReport("Mapping Quality", qc.mappingQuality);
		}
		if (qc.baseQuality != null) {
			logReport("Base Quality", qc.baseQuality);
		}
		if (qc.baq != null) {
			logReport("BAQ", qc.baq);
		}
	}

	private void logReport(String subject, IntegerCounterSet sum) {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		pw.write(subject + " ");
		sum.printReport(pw);
		pw.close();
		String report = sw.toString();
		String[] reportLines = report.split("\\n");
		for (String line : reportLines)
			logger.info(line);
	}


}
