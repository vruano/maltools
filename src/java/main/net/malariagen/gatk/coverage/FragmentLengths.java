package net.malariagen.gatk.coverage;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public abstract class FragmentLengths {

	public class FirstMate {
		int start;
		int stop;
		int length;
		int trim;
	}
	
	private long size;

	protected int maxLength;
	protected Collection<? extends Sample> samples;
	protected Collection<? extends SAMReadGroupRecord> readGroups;
	protected Map<String,Integer> rgIndex;
	protected Map<String,Integer> smIndex;
	protected Map<String,FirstMate> firstMates = new HashMap<String,FirstMate>(100);

	public static long FREQ_SIZE_THR = 1000000;
	
	public static FragmentLengths create(Collection<? extends Sample> samples, Collection<? extends SAMReadGroupRecord> rgs, int maxLength) {
		return new FragmentLengthArrays(samples,rgs,maxLength);
	}
	
	public static FragmentLengths merge(FragmentLengths fl1, FragmentLengths fl2) {
		long newSize = fl1.size + fl2.size;
		
		if (newSize > FREQ_SIZE_THR && !(fl1 instanceof FragmentLengthFrequencies)) {
			FragmentLengthFrequencies result = new FragmentLengthFrequencies(fl1.samples,fl1.readGroups,fl1.maxLength);
			result.mergeIn(fl1);
			result.mergeIn(fl2);
			return result;
		}
		else if (fl1.getClass().equals(fl2.getClass())) {
			fl1.mergeIn(fl2);
			return fl1;
		}
		else if (fl1 instanceof FragmentLengthFrequencies) {
			fl1.mergeIn(fl2);
			return fl1;
		}
		else {
			FragmentLengthFrequencies result = new FragmentLengthFrequencies(fl1.samples,fl1.readGroups,fl1.maxLength);
			result.mergeIn(fl1);
			result.mergeIn(fl2);
			return result;
		}
	}

	
	protected FragmentLengths(Collection<? extends Sample> samples, 
			Collection<? extends SAMReadGroupRecord> rgs, int maxLength) {
		if (maxLength <= 0)
			throw new IllegalArgumentException("max-fragment length must be more than 0");
		this.maxLength = maxLength;
		this.samples = samples;
		this.readGroups = rgs;
		int nextSmIdx = 0;
		int nextRgIdx = 0;
		smIndex = new HashMap<String,Integer>(samples.size());
		rgIndex = new HashMap<String,Integer>(rgs.size());
		for (Sample s : samples) 
			smIndex.put(s.getID(),nextSmIdx++);
		for (SAMReadGroupRecord rg : rgs) 
			smIndex.put(rg.getId(),nextRgIdx++);
		firstMates = new HashMap<String,FirstMate>(100);
		
	}
	
	
	public boolean add(SAMRecord read) {
		if (read.getFirstOfPairFlag())
			return countFirst(read);
		else 
			return countSecond(read);
	}
	
	private boolean countFirst(SAMRecord read) {
		if (!read.getMateUnmappedFlag())
			return false;
		String mateRef = read.getMateReferenceName();
		if (! mateRef.equals("=") && !mateRef.equals(read.getReferenceName()) 
				&& !read.getMateReferenceName().equals("="))
			return false;
		
		int mateStart = read.getMateAlignmentStart();
		int start = read.getAlignmentStart();
		if (start == 0) return false;
		if ((mateStart -start) > maxLength) return false;
		
		FirstMate fm = new FirstMate();
		fm.length = read.getReadLength();
		fm.start = start;
		fm.stop = read.getAlignmentEnd();
		fm.trim =  (fm.stop - fm.start + 1 == fm.length) ? 0 : calculateTrim(read);
		
		firstMates.put(read.getReadName(), fm);
		return false;
	}
	
	private int calculateTrim(SAMRecord read) {
		List<AlignmentBlock> blocks = read.getAlignmentBlocks();
		AlignmentBlock lastBlock = blocks.get(blocks.size() -1);
		return read.getReadLength() - (lastBlock.getReadStart() + lastBlock.getLength());
	}
	
	protected boolean countSecond(SAMRecord read) {
		FirstMate fm = firstMates.remove(read.getReadName());
	    if (fm == null) return false;
	    int start = read.getAlignmentStart();
	    int length = read.getReadLength();
	    int fragmentLength = start - fm.start + 1;
	    int insertLength = fragmentLength - length - fm.length;
	    SAMReadGroupRecord rg = read.getReadGroup();
	    Integer smIdx = smIndex.get(rg.getSample());
	    Integer rgIdx = rgIndex.get(rg.getId());
	    if (smIdx == null)
	    	return false;
	    if (rgIdx == null)
	    	return false;
	    addLengths(fragmentLength, insertLength, smIdx, rgIdx);
	    size++;
	    return true;
	}


	protected abstract void addLengths(int fragmentLength, int insertLength,
			Integer smIdx, Integer rgIdx);

	public void mergeIn(FragmentLengths other) {
		if (other == null)
			return;
		else if (other instanceof FragmentLengthFrequencies) {
			mergeIn((FragmentLengthFrequencies) other);
		}
		else if (other instanceof FragmentLengthArrays) {
			mergeIn((FragmentLengthArrays) other);
		}
		else {
			throw new UnsupportedOperationException();
		}
	}
	
	public void mergeIn(FragmentLengthArrays other) {
		throw new UnsupportedOperationException("yet not implemented");
	}
	
	public void mergeIn(FragmentLengthFrequencies other) {
		throw new UnsupportedOperationException("yet not implemented");
	}

	public static FragmentLengths add(GATKSAMRecord value, FragmentLengths sum) {
		if (sum.add(value)) 
			if (sum instanceof FragmentLengthArrays  && sum.size > FREQ_SIZE_THR) {
				FragmentLengthFrequencies result = new FragmentLengthFrequencies(sum.samples,sum.readGroups,sum.maxLength);
				result.mergeIn(sum);
				return result;
			}
		return sum;
		
	}
}


