package net.malariagen.gatk.coverage;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

public class FragmentLengths {

	public class FirstMate {
		int start;
		int stop;
		int length;
		int trim;
	}

	private int maxLength = 10000;
	
	private Map<String,Integer> rgIndex;
	private Map<String,Integer> smIndex;
	private Map<String,FirstMate> firstMates = new HashMap<String,FirstMate>(100);
	private PriorityQueue<Integer> insertLengths = new PriorityQueue<Integer>(); 
	
	
	private long[][] rgFlFreq;
	private long[][] smFlFreq;
	
	private long[][] rgIlFreq;
	private long[][] smIlFreq;	

	private PriorityQueue<Integer> fragmentLengths = new PriorityQueue<Integer>();
	
	public FragmentLengths(Collection<? extends Sample> samples, Collection<? extends SAMReadGroupRecord> rgs, int maxLength) {
		if (maxLength <= 0)
			throw new IllegalArgumentException("max-fragment length must be more than 0");
		this.maxLength = maxLength;
		int nextSmIdx = 0;
		int nextRgIdx = 0;
		smIndex = new HashMap<String,Integer>(samples.size());
		rgIndex = new HashMap<String,Integer>(rgs.size());
		for (Sample s : samples) 
			smIndex.put(s.getID(),nextSmIdx++);
		for (SAMReadGroupRecord rg : rgs) 
			smIndex.put(rg.getId(),nextRgIdx++);
		
		rgFlFreq = new long[nextRgIdx][maxLength];
		smFlFreq = new long[nextSmIdx][maxLength];
		rgIlFreq = new long[nextRgIdx][maxLength];
		smIlFreq = new long[nextSmIdx][maxLength];
	}
	
	public void count(SAMRecord read) {
		if (read.getFirstOfPairFlag())
			countFirst(read);
		else 
			countSecond(read);
	}
	
	public void countFirst(SAMRecord read) {
		if (!read.getMateUnmappedFlag())
			return;
		String mateRef = read.getMateReferenceName();
		if (! mateRef.equals("=") && !mateRef.equals(read.getReferenceName()) 
				&& !read.getMateReferenceName().equals("="))
			return;
		
		int mateStart = read.getMateAlignmentStart();
		int start = read.getAlignmentStart();
		if (start == 0) return;
		if ((mateStart -start) > maxLength) return;
		
		FirstMate fm = new FirstMate();
		fm.length = read.getReadLength();
		fm.start = start;
		fm.stop = read.getAlignmentEnd();
		fm.trim =  (fm.stop - fm.start + 1 == fm.length) ? 0 : calculateTrim(read);
	}
	
	public int calculateTrim(SAMRecord read) {
		List<AlignmentBlock> blocks = read.getAlignmentBlocks();
		AlignmentBlock lastBlock = blocks.get(blocks.size() -1);
		return read.getReadLength() - (lastBlock.getReadStart() + lastBlock.getLength());
	}
	
	public void countSecond(SAMRecord read) {
		FirstMate fm = firstMates.remove(read.getReadName());
	    if (fm == null) return;
	    int start = read.getAlignmentStart();
	    int length = read.getReadLength();
	    int fragmentLength = start - fm.start + 1;
	    int insertLength = fragmentLength - length - fm.length;
	    SAMReadGroupRecord rg = read.getReadGroup();
	    Integer smIdx = smIndex.get(rg.getSample());
	    Integer rgIdx = rgIndex.get(rg.getId());
	    if (smIdx == null)
	    	return;
	    if (rgIdx == null)
	    	return;
	    smIlFreq[smIdx][0]++;
	    rgIlFreq[rgIdx][0]++;
	    smFlFreq[smIdx][0]++;
	    rgFlFreq[rgIdx][0]++;
	    
	    if (smIlFreq[smIdx][insertLength]++ == 0) 
	    	insertLengths.offer(insertLength);
	    if (smFlFreq[smIdx][fragmentLength]++ == 0)
	    	fragmentLengths.offer(fragmentLength);
	    rgIlFreq[rgIdx][insertLength]++;
	    rgFlFreq[rgIdx][fragmentLength]++;
	}
	
	public void mergeIn(FragmentLengths other) {
		if (other == null)
			return;
		if (other.maxLength > this.maxLength) 
			increaseMaxLength(other.maxLength);
		int size = other.insertLengths.size();
		
		if (size << 1 > maxLength) 
			arrayBasedIlMerge(other);
		else
			pqBasedIlMerge(other);
	}

	private void pqBasedIlMerge(FragmentLengths other) {
		// TODO Auto-generated method stub
		
	}

	private void arrayBasedIlMerge(FragmentLengths other) {
		for (Map.Entry<String,Integer> e : smIndex) {
			String s = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.smIndex.get(s);
			if (otIdxI == null)
				continue;
			arrayBasedIlMerge(other, myIdx, otIdxI);
		}
		for (Map.Entry<String,Integer> e: rgIndex) {
			String rg = e.getKey();
		}
	}

	private void arrayBasedIlMerge(FragmentLengths other, int myIdx,
			Integer otIdxI) {
		
		// TODO Auto-generated method stub
		
	}

	private void increaseMaxLength(int newValue) {
		for (int i = 0; i < rgIlFreq.length; i++) {
			rgIlFreq[i] = Arrays.copyOf(rgIlFreq[i],newValue);
			rgFlFreq[i] = Arrays.copyOf(rgFlFreq[i],newValue);
		}
		for (int i = 0; i < smIlFreq.length; i++) {
			smIlFreq[i] = Arrays.copyOf(smIlFreq[i],newValue);
			smIlFreq[i] = Arrays.copyOf(smIlFreq[i],newValue);
		}
	}
	
	
}