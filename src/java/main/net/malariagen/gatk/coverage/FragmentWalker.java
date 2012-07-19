package net.malariagen.gatk.coverage;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@ReadFilters({ BadMateFilter.class, UnmappedReadFilter.class,
	NotPrimaryAlignmentFilter.class })
@By(DataSource.READS)
@Requires(DataSource.READS)
@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
public abstract class FragmentWalker<MapType,ReduceType> extends ReadWalker<GATKSAMRecord,FragmentWalkerReduceType<ReduceType>> {

	@Argument(shortName="maxfl", fullName="max_fragment_length", 
			  doc = "maximum admissible fragment length; mapped read pairs on the same chromosome that have a larger fragment size would be discarded. Large max-fragment length results in a increased need of memory to keep record of the first read until its mate is found", required = false)
	public int maxLength = 1000;
	
	@Argument(shortName="mmq", fullName="min_mapping_quality",
			doc = "minimum mapping quality", required = false)
	public int minimumMappingQuality = 0;

	
	private IndexedFastaSequenceFile reference;

	private GenomeLocParser locParser;

	private SAMSequenceDictionary dictionary;
		
	public abstract ReduceType reduceFragmentInit();
	public abstract MapType mapFragment(ReferenceContext ref, FragmentRecord fragment, FragmentMetaDataTracker metaDataTracker);
	public abstract ReduceType reduceFragment(MapType value,ReduceType sum);
	
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
	}


	@Override
	public FragmentWalkerReduceType<ReduceType> reduceInit() {
		return new FragmentWalkerReduceType<ReduceType>(maxLength, minimumMappingQuality,reduceFragmentInit());
	}

	@Override
	public FragmentWalkerReduceType<ReduceType> reduce(
			GATKSAMRecord value,
			FragmentWalkerReduceType<ReduceType> sum) {
		FragmentRecord frecord = sum.add(value);
		if (frecord != null) {
			GenomeLoc loc = frecord.getRegion(locParser);
			ReferenceContext ref = new ReferenceContext(locParser, loc, loc, buildBases(loc));
			FragmentMetaDataTracker tracker = new FragmentMetaDataTracker();
			MapType mtValue = this.mapFragment(ref, frecord, tracker);
			sum.setReduceObject(reduceFragment(mtValue,sum.getReduceObject()));
		}
		return sum;
	}
	
	private byte[] buildBases(GenomeLoc loc) {
        
        long start = loc.getStart();
        long stop = loc.getStop();

        // Read with no aligned bases?  Return an empty array.
        if(stop - start + 1 == 0)
            return new byte[0];

    
        try {
          SAMSequenceRecord srecord = dictionary.getSequence(loc.getContig());
          if (srecord.getSequenceLength() < stop)
        	  stop = srecord.getSequenceLength();
          ReferenceSequence subsequence = reference.getSubsequenceAt(loc.getContig(), start, stop);    
          return subsequence.getBases();
        }
        catch (RuntimeException e) {
        	System.err.println(loc);
        	throw e;
        }
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
		this.onFragmentTravesalDone(result.getReduceObject());
	}

	@Override
	public boolean isReduceByInterval() {
		return false;
	}
	
	public void onFragmentTravesalDone(ReduceType sum) {
		
	}
	
	

}
