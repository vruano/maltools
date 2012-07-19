package net.malariagen.gatk.coverage;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;


public class FragmentLengthArrays extends FragmentLengths {

	int[][] smIlengths;
	int[][] rgIlengths;
	int[][] smFlengths;
	int[][] rgFlengths;
	
	
	FragmentLengthArrays(Collection<String> samples,
			Collection<String> rgs, int maxLength, int minimumMappingQuality) {
		super(samples, rgs, maxLength, minimumMappingQuality);
		int capacity = 10000;
		int sampleCount = smIndex.size();
		int rgCount = rgIndex.size();
		smIlengths = new int[sampleCount][capacity + 1];
		smFlengths = new int[sampleCount][capacity + 1];		
		rgIlengths = new int[rgCount][capacity + 1];
		rgFlengths = new int[rgCount][capacity + 1];	
	}

	

	@Override
	protected void addLengths(int fragmentLength, int insertLength,
			Integer smIdxI, Integer rgIdxI) {
		int smIdx = smIdxI;
		int rgIdx = rgIdxI;
		int smNext;
		int rgNext;
		if ((smNext = ++smIlengths[smIdx][0]) >= smIlengths[smIdx].length) {
			smIlengths[smIdx] = Arrays.copyOf(smIlengths[smIdx], smIlengths[smIdx].length << 1);
			smFlengths[smIdx] = Arrays.copyOf(smFlengths[smIdx], smFlengths[smIdx].length << 1);
		}
		if ((rgNext = ++rgIlengths[rgIdx][0]) >= rgIlengths[rgIdx].length) {
			rgIlengths[rgIdx] = Arrays.copyOf(rgIlengths[rgIdx], rgIlengths[rgIdx].length << 1);
			rgFlengths[rgIdx] = Arrays.copyOf(rgFlengths[rgIdx], rgFlengths[rgIdx].length << 1);
		}
		smFlengths[smIdx][0]++;
		rgFlengths[rgIdx][0]++;
		smIlengths[smIdx][smNext] = rgIlengths[rgIdx][rgNext] = insertLength;
		smFlengths[smIdx][smNext] = rgFlengths[rgIdx][rgNext] = fragmentLength;
	}

	@Override
	public void mergeIn(FragmentLengthArrays other) {
		mergeSampleArrays(other);
		mergeReadGroupArrays(other);
		this.size += other.size;
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
			int otCount = otArray[0];
			int myCount = this.rgIlengths[myIdx][0];
			if (otCount + myCount + 1 >= rgIlengths[myIdx].length) 
				rgIlengths[myIdx] = Arrays.copyOf(rgIlengths[myIdx],(otCount + myCount + 1) * 2);
			System.arraycopy(otArray,1,rgIlengths[myIdx],myCount + 1,otCount);
			rgIlengths[myIdx][0] += otCount;
			
			otArray = other.rgFlengths[otIdx];
			otCount = otArray[0];
			myCount = this.rgFlengths[myIdx][0];
			if (otCount + myCount + 1 >= rgFlengths[myIdx].length) 
				rgFlengths[myIdx] = Arrays.copyOf(rgFlengths[myIdx],(otCount + myCount + 1) * 2);
			System.arraycopy(otArray,1,rgFlengths[myIdx],myCount + 1,otCount);
			rgFlengths[myIdx][0] += otCount;
		}
	}

	private void mergeSampleArrays(FragmentLengthArrays other) {
		for (Map.Entry<String,Integer> e: smIndex.entrySet()) {
			String rg = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.smIndex.get(rg);
			if (otIdxI == null)
				continue;
			int otIdx = otIdxI;
			int[] otArray = other.smIlengths[otIdx];
			int otCount = otArray[0];
			int myCount = this.smIlengths[myIdx][0];
			if (otCount + myCount + 1 >= smIlengths[myIdx].length) 
				smIlengths[myIdx] = Arrays.copyOf(smIlengths[myIdx],(otCount + myCount + 1) * 2);
			System.arraycopy(otArray,1,smIlengths[myIdx],myCount + 1,otCount);
			smIlengths[myIdx][0] += otCount;
			
			otArray = other.smFlengths[otIdx];
			otCount = otArray[0];
			myCount = this.smFlengths[myIdx][0];
			if (otCount + myCount + 1 >= smFlengths[myIdx].length) 
				smFlengths[myIdx] = Arrays.copyOf(smFlengths[myIdx],(otCount + myCount + 1) * 2);
			System.arraycopy(otArray,1,smFlengths[myIdx],myCount + 1,otCount);
			smFlengths[myIdx][0] += otCount;
		}
	}

	public FragmentLengthFrequencies frequencies() {
		FragmentLengthFrequencies result = new FragmentLengthFrequencies(this.samples,this.readGroups,this.maxLength,this.minimumMappingQuality);
		result.mergeIn(this);
		result.listeners.addAll(this.listeners);
		return result;
	}

	@Override
	public FragmentLengthSummary summary() {
		return frequencies().summary();
	}

}
