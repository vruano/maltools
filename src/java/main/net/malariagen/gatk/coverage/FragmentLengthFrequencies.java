package net.malariagen.gatk.coverage;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Stack;

import org.broadinstitute.sting.gatk.samples.Sample;
import net.sf.samtools.SAMReadGroupRecord;

public class FragmentLengthFrequencies extends FragmentLengths {

	
	private PriorityQueue<Integer> insertLengths = new PriorityQueue<Integer>(); 
	
	
	private long[][] rgFlFreq;
	private long[][] smFlFreq;
	
	private long[][] rgIlFreq;
	private long[][] smIlFreq;	

	private PriorityQueue<Integer> fragmentLengths = new PriorityQueue<Integer>();

	
	
	FragmentLengthFrequencies(Collection<? extends Sample> samples, Collection<? extends SAMReadGroupRecord> rgs, int maxLength) {
		super(samples,rgs,maxLength);
		int sampleCount = smIndex.size();
		int rgCount = rgIndex.size();
		rgFlFreq = new long[rgCount][maxLength];
		smFlFreq = new long[sampleCount][maxLength];
		rgIlFreq = new long[rgCount][maxLength];
		smIlFreq = new long[sampleCount][maxLength];
	}


	@Override
	protected void addLengths(int fragmentLength, int insertLength,
			Integer smIdx, Integer rgIdx) {
		smIlFreq[smIdx][0]++;
	    rgIlFreq[rgIdx][0]++;
	    smFlFreq[smIdx][0]++;
	    rgFlFreq[rgIdx][0]++;
	    
	    if (smIlFreq[smIdx][insertLength]++ == 0) 
	    	insertLengths.add(insertLength);
	    if (smFlFreq[smIdx][fragmentLength]++ == 0)
	    	fragmentLengths.add(fragmentLength);
	    rgIlFreq[rgIdx][insertLength]++;
	    rgFlFreq[rgIdx][fragmentLength]++;
	}
	

	
	public void mergeIn(FragmentLengthArrays other) {
		if (other.maxLength > this.maxLength)
			increaseMaxLength(other.maxLength);
		mergeSampleArrays(other);
		mergeReadGroupArrays(other);
	}
	
	private void mergeReadGroupArrays(FragmentLengthArrays other) {
		for (Map.Entry<String,Integer> e: rgIndex.entrySet()) {
			String rg = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.rgIndex.get(rg);
			if (otIdxI == null)
				continue;
			int otIdx = otIdxI;
			int[] otArray = other.rgIlengths[otIdx];
			int count = otArray[0];
			long[] myFreq = this.rgIlFreq[myIdx];
			for (int i = 1; i <= count; i++)
				myFreq[otArray[i]]++;
			myFreq[0] += count;
			otArray = other.rgFlengths[otIdx];
			count = otArray[0];
			myFreq = this.rgFlFreq[myIdx];
			for (int i = 1; i <= count; i++)
				myFreq[otArray[i]]++;
			myFreq[0] += count;			
		}
	}


	private void mergeSampleArrays(FragmentLengthArrays other) {
		for (Map.Entry<String,Integer> e : smIndex.entrySet()) {
			String s = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.smIndex.get(s);
			if (otIdxI == null)
				continue;
			int otIdx = otIdxI;
			int[] otArray = other.smIlengths[otIdx];
			int count = otArray[0];
			long[] myFreq = this.smIlFreq[myIdx];
			for (int i = 1; i <= count; i++)
				myFreq[otArray[i]]++;
			myFreq[0] += count;
			otArray = other.smFlengths[otIdx];
			count = otArray[0];
			myFreq = this.smFlFreq[myIdx];
			for (int i = 1; i <= count; i++)
				myFreq[otArray[i]]++;
			myFreq[0] += count;			
		}
	}


	public void mergeIn(FragmentLengthFrequencies other) {
		if (other == null)
			return;
		if (other.maxLength > this.maxLength) 
			increaseMaxLength(other.maxLength);
		int size = other.insertLengths.size();
			mergeFreqs(other, size << 1 > maxLength);
	}

	private void mergeFreqs(FragmentLengthFrequencies other, boolean arrayBased) {
		mergeSampleFreqs(other, arrayBased);
		mergeReadGroupFreqs(other, arrayBased);
	}

	private void mergeReadGroupFreqs(FragmentLengthFrequencies other, boolean arrayBased) {
		for (Map.Entry<String,Integer> e: rgIndex.entrySet()) {
			String rg = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.rgIndex.get(rg);
			if (otIdxI == null)
				continue;
			if (arrayBased) {
				arrayBasedMerge(other, myIdx, otIdxI,true,true);
				arrayBasedMerge(other, myIdx, otIdxI,true,false);			
			}
			else {
				pgBasedMerge(other, myIdx, otIdxI,true,true);
				pgBasedMerge(other, myIdx, otIdxI,true,false);
			}
		}
	}

	private void mergeSampleFreqs(FragmentLengthFrequencies other, boolean arrayBased) {
		for (Map.Entry<String,Integer> e : smIndex.entrySet()) {
			String s = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.smIndex.get(s);
			if (otIdxI == null)
				continue;
			if (arrayBased) {
			arrayBasedMerge(other, myIdx, otIdxI,false,true);
			arrayBasedMerge(other, myIdx, otIdxI,false,false);			
			}
			else {
				pgBasedMerge(other, myIdx, otIdxI,false,true);
				pgBasedMerge(other, myIdx, otIdxI,false,false);
			}
		}
	}

	private void pgBasedMerge(FragmentLengthFrequencies ot, int myIdx, int otIdx,
			boolean rg, boolean il) {
		PriorityQueue<Integer> pq = il ? ot.insertLengths : ot.fragmentLengths;
		Stack<Integer> st = new Stack<Integer>();
		st.ensureCapacity(pq.size());
		long[] myArray = rg ? (il ? rgIlFreq[myIdx] : rgFlFreq[myIdx]) : (il ? smIlFreq[myIdx] : smFlFreq[myIdx]);
		long[] otArray = rg ? (il ? ot.rgIlFreq[otIdx] : ot.rgFlFreq[otIdx]) : (il ? ot.smIlFreq[otIdx] : ot.smFlFreq[otIdx]);
		while (!pq.isEmpty()) {
			Integer next = pq.poll();
			myArray[next] += otArray[next];
			st.push(next);
		}
		while (!st.isEmpty())
			pq.add(st.pop());
		myArray[0] += otArray[0];
	}

	private void arrayBasedMerge(FragmentLengthFrequencies ot, int myIdx,
			int otIdx, boolean rg, boolean il) {
		long[] myArray = rg ? (il ? rgIlFreq[myIdx] : rgFlFreq[myIdx]) : (il ? smIlFreq[myIdx] : smFlFreq[myIdx]);
		long[] otArray = rg ? (il ? ot.rgIlFreq[otIdx] : ot.rgFlFreq[otIdx]) : (il ? ot.smIlFreq[otIdx] : ot.smFlFreq[otIdx]);
		for (int i = 0; i < otArray.length; i++)
			myArray[i] += otArray[i];
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