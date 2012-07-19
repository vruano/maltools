package net.malariagen.gatk.coverage;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class FragmentWalkerReduceType<ReduceType> {

	class FirstMate {
		GATKSAMRecord record;
		int length;
		int start;
		int stop;
		int trim;
		
		FirstMate(GATKSAMRecord read) {
			record = read;
			length = record.getReadLength();
			start = record.getAlignmentStart();
			stop = record.getAlignmentEnd();
			trim = (stop - start + 1 == length) ? 0 :  calculateTrim(record);
		}
	
	}

	protected long size;

	protected int maxLength;
	protected Map<String, FirstMate> firstMates;

	protected int minimumMappingQuality = 0;
	
	protected ReduceType sum;

	FragmentWalkerReduceType(int maxLength, int minimumMappingQuality, ReduceType sum) {
		this.maxLength = maxLength;
		this.minimumMappingQuality = minimumMappingQuality;
		firstMates = new HashMap<String, FirstMate>(100);
		this.sum = sum;

	}
	
	ReduceType getReduceObject() {
		return sum;
	}
	
	void setReduceObject(ReduceType o) {
		sum = o;
	}
	

	
	void mergeIn(FragmentWalkerReduceType<ReduceType> other) {
		this.firstMates.putAll(other.firstMates);
	}
	
	public FragmentRecord add(GATKSAMRecord read) {
		if (read == null)
			throw new IllegalArgumentException("input read cannot be null");
		FirstMate fm =  firstMates.remove(read.getReadName());;
		if (fm != null)
			return countSecond(read,fm);
		else {
			countFirst(read);
			return null;
		}
	}

	private void countFirst(GATKSAMRecord read) {
		if (read.getReadUnmappedFlag()) return;
		if (read.getMateUnmappedFlag()) return;
		if (read.getMappingQuality() < minimumMappingQuality) return;
		String mateRef = read.getMateReferenceName();
		if (!mateRef.equals("=") && !mateRef.equals(read.getReferenceName())) return;
		int mateStart = read.getMateAlignmentStart();
		int start = read.getAlignmentStart();
		if (mateStart - start > maxLength || start - mateStart > maxLength) return;
		if (read.getReadNegativeStrandFlag()) return;
		FirstMate fm = new FirstMate(read);
		firstMates.put(read.getReadName(), fm);
	}

	private int calculateTrim(SAMRecord read) {
		List<AlignmentBlock> blocks = read.getAlignmentBlocks();
		int lastPos = 0;
		int result = 0;
		for (AlignmentBlock b : blocks) {
			result += b.getReadStart() - lastPos - 1;
			lastPos = b.getReadStart() + b.getLength() - 1;
		}
		result += read.getReadLength() - lastPos - 1;
		return result;
	}

	protected FragmentRecord countSecond(GATKSAMRecord read, FirstMate fm) {
		if (!read.getReadNegativeStrandFlag())
			return null;
		if (read.getMappingQuality() < minimumMappingQuality)
			return null;
		int start = read.getAlignmentStart();
		int end = read.getAlignmentEnd();
		int length = read.getReadLength();
		int trim = ((end - start + 1 == length) ? 0 : calculateTrim(read));
		int fragmentLength = end - fm.start + 1 + fm.trim + trim;
		SAMReadGroupRecord rg = read.getReadGroup();
		if (rg == null)
			throw new IllegalArgumentException(
					"reads need to belong to a read-group");
		if (fragmentLength <= fm.length + length)
			return null;
		if (fragmentLength > maxLength)
			return null;
		
		FragmentRecord result;
		if (read.getReadNegativeStrandFlag()) 
			result = new FragmentRecord(read,fm.record,fragmentLength);
		else
			result = new FragmentRecord(fm.record,read,fragmentLength);
		return result;
	}

	public void mergeIn(FragmentLengthArrays other) {
		throw new UnsupportedOperationException("yet not implemented");
	}

	public void mergeIn(FragmentLengthFrequencies other) {
		throw new UnsupportedOperationException("yet not implemented");
	}

	public long size() {
		return size;
	}


}
