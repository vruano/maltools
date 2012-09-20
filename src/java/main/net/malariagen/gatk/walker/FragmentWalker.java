package net.malariagen.gatk.walker;

import java.util.HashSet;

import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import net.malariagen.gatk.filters.UnmappedMateFilter;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@ReadFilters({ NotPrimaryAlignmentFilter.class })
@By(DataSource.READS)
@Requires(DataSource.READS)
@Downsample(by=DownsampleType.NONE)
@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
public abstract class FragmentWalker<MapType, ReduceType> extends
		ReadWalker<GATKSAMRecord, FragmentWalkerReduceType<ReduceType>> {

	@Argument(shortName = "maxfl", fullName = "maximumFragmentLength", doc = "maximum admissible fragment length; mapped read pairs on the same chromosome that have a larger fragment size would be discarded. Large max-fragment length results in a increased need of memory to keep record of the first read until its mate is found", required = false)
	public int maxLength = 10000;

	@Argument(shortName = "mmq", fullName = "minimumMappingQuality", doc = "minimum mapping quality", required = false)
	public int minimumMappingQuality = 0;

	@Argument(shortName = "filter", doc = "filters to apply", required = false)
	protected Set<FragmentFilter> filters = null;

	@Argument(shortName = "sort", doc = "how to sort the fragments", required = false)
	protected FragmentOrder order;

	private SortedMap<GenomeLoc, List<FragmentRecord>> sortedBuffer;
	private SortedMap<GenomeLoc, Integer> pendingStarts;

	private IndexedFastaSequenceFile reference;

	GenomeLocParser locParser;

	private SAMSequenceDictionary dictionary;

	public abstract ReduceType reduceFragmentInit();

	public abstract MapType mapFragment(ReferenceContext ref,
			FragmentRecord fragment, FragmentMetaDataTracker metaDataTracker);

	public abstract ReduceType reduceFragment(MapType value, ReduceType sum);

	@Override
	public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read,
			ReadMetaDataTracker metaDataTracker) {
		return read;
	}

	@Override
	public void initialize() {
		super.initialize();
		reference = getToolkit().getReferenceDataSource().getReference();
		dictionary = getToolkit().getMasterSequenceDictionary();
		locParser = getToolkit().getGenomeLocParser();
		Set<FragmentFilter> finalFilters = filters == null ? new HashSet<FragmentFilter>() : filters;
		if (filters == null) {
			filters = new HashSet<FragmentFilter>();
			filters.add(FragmentFilter.UNMAPPED);
			filters.add(FragmentFilter.TRANSLOCATION);
			filters.add(FragmentFilter.INVERSION);
		}
		if (minimumMappingQuality > 0 && filters == null) {
			filters.add(FragmentFilter.UNMAPPED);
			filters.add(FragmentFilter.LOW_MQ);
		}
		filters = finalFilters;
		if (order == null)
			order = FragmentOrder.NONE;
		if (order != FragmentOrder.NONE) {
			sortedBuffer = new TreeMap<GenomeLoc, List<FragmentRecord>>();
			pendingStarts = new TreeMap<GenomeLoc,Integer>();
		}
	}

	@Override
	public FragmentWalkerReduceType<ReduceType> reduceInit() {
		return new FragmentWalkerReduceType<ReduceType>(this,
				reduceFragmentInit());
	}

	@Override
	public FragmentWalkerReduceType<ReduceType> reduce(GATKSAMRecord value,
			FragmentWalkerReduceType<ReduceType> sum) {
		if (!passFilters(value))
			return sum;
		FragmentRecord frecord = sum.add(value);
		if (frecord == null) {
			if (order != FragmentOrder.NONE) addToPendingStarts(value);
		} else {
			if (order == FragmentOrder.NONE) {
				if (passFilters(frecord))
					processAndReduceFragment(sum, frecord);
			} else {
				if (!addToSortedBuffer(frecord)) 
					processAndReduceFragment(sum, frecord);
				else 
				    processAndReduceSortedFragments(sum);
			}
		}
		return sum;
	}


	private boolean passFilters(GATKSAMRecord value) {
		if (!value.getReadPairedFlag())
			return false;
		return true;
//		if (this.maxLength <= 0)
//			return true;
//		if (value.getReadUnmappedFlag())
//			return true;
//		int myStart = value.getAlignmentStart();
//		int otherStart = value.getMateAlignmentStart();
//		int firstStart = Math.min(myStart, otherStart);
//		int secondStart = Math.max(myStart,otherStart);
//		return (secondStart - firstStart + value.getReadLength() - 1 <= maxLength); 
	}


	private boolean passFilters(FragmentRecord result) {

		for (FragmentFilter filter : this.filters) {
			switch (filter) {
			case UNMAPPED:
				if (result.getFirstRead().getReadUnmappedFlag())
					return false;
				if (result.getLastRead().getReadUnmappedFlag())
					return false;
				break;
			case LAST_UNMAPPED:
				if (result.getLastRead().getReadUnmappedFlag())
					return false;
				break;
			case FIRST_UNMAPPED:
				if (result.getFirstRead().getReadUnmappedFlag())
					return false;
				break;
			case BOTH_UNMAPPED:
				if (result.getFirstRead().getReadUnmappedFlag()
						&& result.getLastRead().getReadUnmappedFlag())
					return false;
				break;
			case TRANSLOCATION:
				if (result.isTranslocation(this.maxLength))
					return false;
				break;
			case INVERSION:
				if (!result.isMapped()
						|| result.isTranslocation(this.maxLength))
					break;
				if (result.getReferenceFirstRead().getReadNegativeStrandFlag()
						|| !result.getReferenceLastRead()
								.getReadNegativeStrandFlag())
					return false;
				break;
			case BOTH_FACE_AWAY:
				if (!result.isMapped()
						|| result.isTranslocation(this.maxLength))
					break;
				if (result.getReferenceFirstRead().getReadNegativeStrandFlag()
						&& !result.getReferenceLastRead()
								.getReadNegativeStrandFlag())
					return false;
				break;
			case FIRST_FACE_AWAY:
				if (!result.isMapped()
						|| result.isTranslocation(this.maxLength))
					break;
				if (result.getReferenceFirstRead().getReadNegativeStrandFlag())
					return false;
				break;
			case LAST_FACE_AWAY:
				if (!result.isMapped()
						|| result.isTranslocation(this.maxLength))
					break;
				if (!result.getReferenceLastRead().getReadNegativeStrandFlag())
					return false;
				break;
			case LOW_MQ:
				if (result.getFirstRead().getMappingQuality() < this.minimumMappingQuality)
					return false;
				if (result.getLastRead().getMappingQuality() < this.minimumMappingQuality)
					return false;
				break;
			case BOTH_LOW_MQ:
				if (result.getFirstRead().getMappingQuality() < this.minimumMappingQuality
						&& result.getLastRead().getMappingQuality() < this.minimumMappingQuality)
					return false;
				break;
			case FIRST_LOW_MQ:
				if (result.getFirstRead().getMappingQuality() < this.minimumMappingQuality)
					return false;
				break;
			case LAST_LOW_MQ:
				if (result.getLastRead().getMappingQuality() < this.minimumMappingQuality)
					return false;
				break;
			default:
			}
		}
		return true;
	}

	private boolean sortable(GATKSAMRecord value) {
		if (value.getReadUnmappedFlag()) return false;
		if (value.getMateUnmappedFlag()) return false;
		if (value.getReferenceIndex() != value.getMateReferenceIndex()) return false;
		if (order != FragmentOrder.REFERENCE_START)
		  if (!value.getFirstOfPairFlag()) return false;
		if (this.maxLength < value.getMateAlignmentStart() + value.getReadLength() - value.getAlignmentStart() +1)
		  return false;
		return true;
	}

	private boolean addToPendingStarts(GATKSAMRecord value) {
		GenomeLoc loc = startLocationFrom(value);
		if (loc == null)
			return false;
		else {
			addToPendingStarts(loc);
			return true;
		}
	}
	
   private void addToPendingStarts(GenomeLoc loc) {
		
		Integer count = pendingStarts.get(loc);
		if (count == null)
			pendingStarts.put(loc, 1);
		else
			pendingStarts.put(loc, count + 1);
	}

	private GenomeLoc startLocationFrom(GATKSAMRecord value) {
		if (!sortable(value))
			return null;
		GenomeLoc loc = null;
		if (order == FragmentOrder.FRAGMENT_START) {
			if (value.getReadNegativeStrandFlag())
				loc = locParser.createGenomeLoc(value.getReferenceName(),
						value.getAlignmentEnd(), value.getAlignmentEnd());
			else
				loc = locParser.createGenomeLoc(value.getReferenceName(),value.getAlignmentStart(),value.getAlignmentStart());
		} else
			loc = locParser.createGenomeLoc(value.getReferenceName(),value.getAlignmentStart(),value.getAlignmentStart());
		return loc;
	}
	
	private boolean addToSortedBuffer(FragmentRecord frecord) {
		GenomeLoc pendingLoc = startLocationFrom(frecord.getFirstReadReported());
		GenomeLoc loc;
		if (pendingLoc == null) {
			loc = startLocationFrom(frecord.getLastReadReported());
		}
		else {
			loc = pendingLoc;
		}
		if (loc == null)
			return false;
		List<FragmentRecord> list = sortedBuffer.get(loc);
		if (list == null)
			sortedBuffer.put(loc, list = new LinkedList<FragmentRecord>());
		list.add(frecord);
		if (pendingLoc != null) {
			Integer pendingStartsCount = pendingStarts.get(pendingLoc);
			if (pendingStartsCount == null || pendingStartsCount == 0)
				throw new IllegalStateException("cannot be");
			if (pendingStartsCount == 1)
				pendingStarts.remove(pendingLoc);
			else
				pendingStarts.put(pendingLoc,pendingStartsCount - 1);
		}
		return true;
	}




	private void processAndReduceSortedFragments(
			FragmentWalkerReduceType<ReduceType> sum) {
		GenomeLoc loc = pendingStarts.size() == 0 ? null : pendingStarts
				.firstKey();
		if (loc != null && sortedBuffer.firstKey().compareTo(loc) > 0)
			return;
		SortedMap<GenomeLoc, List<FragmentRecord>> ready = loc == null ? sortedBuffer
				: sortedBuffer.headMap(loc);
		for (List<FragmentRecord> rl : ready.values())
			for (FragmentRecord f : rl) {
				if (passFilters(f)) {
					processAndReduceFragment(sum, f);
				}
			}
		ready.clear();
	}
	
	
	private void processAndReduceFragment(
			FragmentWalkerReduceType<ReduceType> sum, FragmentRecord frecord) {
		SAMRecord referenceFirstRead = frecord.getReferenceFirstRead();
		GenomeLoc loc = frecord.isTranslocation() ? locParser
				.createGenomeLoc(referenceFirstRead.getReferenceName(),
						         referenceFirstRead.getAlignmentStart(),
						         referenceFirstRead.getAlignmentStart()) : frecord
				.getRegion(locParser);
		ReferenceContext ref = new ReferenceContext(locParser, loc, loc,
				buildBases(loc));
		FragmentMetaDataTracker tracker = new FragmentMetaDataTracker();
		MapType mtValue = this.mapFragment(ref, frecord, tracker);
		sum.setReduceObject(reduceFragment(mtValue, sum.getReduceObject()));
	}

	protected boolean requiresStartOrder() {
		return false;
	}

	private byte[] buildBases(GenomeLoc loc) {

		long start = loc.getStart();
		long stop = loc.getStop();

		// Read with no aligned bases? Return an empty array.
		if (stop - start + 1 == 0)
			return new byte[0];

		SAMSequenceRecord srecord = dictionary.getSequence(loc.getContig());
		if (srecord == null)
			return null;
//			throw new IllegalStateException("Invalid loc contig " + loc.getContig());
		if (srecord.getSequenceLength() < stop)
			stop = srecord.getSequenceLength();
		ReferenceSequence subsequence = reference.getSubsequenceAt(
				loc.getContig(), start, stop);
		return subsequence.getBases();
	}

	@Override
	public boolean requiresOrderedReads() {
		return true;
	}

	@Override
	public boolean includeReadsWithDeletionAtLoci() {
		return true;
	}

	@Override
	public void onTraversalDone(FragmentWalkerReduceType<ReduceType> result) {
		super.onTraversalDone(result);
		if (order != FragmentOrder.NONE) {
			if (this.pendingStarts.size() > 0)
				throw new IllegalStateException(
						"inbalance detected some read-pair missing");
			if (this.sortedBuffer.size() > 0)
				throw new IllegalStateException("not all fragments reduced");
		}
		this.onFragmentTravesalDone(result.getReduceObject());
	}

	
	@Override
	public boolean isReduceByInterval() {
		return false;
	}

	public void onFragmentTravesalDone(ReduceType sum) {

	}

}
