package net.malariagen.gatk.walker;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Stack;

public class FragmentLengthFrequencies extends FragmentLengths {

	private PriorityQueue<Integer> insertLengths = new PriorityQueue<Integer>();
	private PriorityQueue<Integer> fragmentLengths = new PriorityQueue<Integer>();
	private boolean[] hasInsertLength;
	private boolean[] hasFragmentLength;

	long[][] rgFlFreq;
	long[][] smFlFreq;

	long[][] rgIlFreq;
	long[][] smIlFreq;

	FragmentLengthFrequencies(Collection<String> samples,
			Collection<String> rgs, int maxLength, int minimumMappingQuality) {
		super(samples, rgs, maxLength, minimumMappingQuality);
		int sampleCount = smIndex.size();
		int rgCount = rgIndex.size();
		hasInsertLength = new boolean[maxLength + 1];
		hasFragmentLength = new boolean[maxLength + 1];
		rgFlFreq = new long[rgCount][maxLength + 1];
		smFlFreq = new long[sampleCount][maxLength + 1];
		rgIlFreq = new long[rgCount][maxLength + 1];
		smIlFreq = new long[sampleCount][maxLength + 1];
	}


	@Override
	protected void addLengths(int fragmentLength, int insertLength,
			Integer smIdx, Integer rgIdx) {
		smIlFreq[smIdx][0]++;
		rgIlFreq[rgIdx][0]++;
		smFlFreq[smIdx][0]++;
		rgFlFreq[rgIdx][0]++;

		if (!hasInsertLength[insertLength]) {
			insertLengths.add(insertLength);
			hasInsertLength[insertLength] = true;
		}
		if (!hasFragmentLength[fragmentLength]) {
			fragmentLengths.add(fragmentLength);
			hasFragmentLength[fragmentLength] = true;
		}
		rgFlFreq[rgIdx][fragmentLength]++;
		rgIlFreq[rgIdx][insertLength]++;
		smIlFreq[smIdx][insertLength]++;
		smFlFreq[smIdx][fragmentLength]++;
	}

	public void mergeIn(FragmentLengthArrays other) {
		if (other.maxLength > this.maxLength)
			increaseMaxLength(other.maxLength);
		mergeSampleArrays(other);
		mergeReadGroupArrays(other);
		this.size += other.size;
	}

	private void mergeReadGroupArrays(FragmentLengthArrays other) {
		for (Map.Entry<String, Integer> e : rgIndex.entrySet()) {
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
				if (myFreq[otArray[i]]++ == 0 && !hasInsertLength[otArray[i]]) {
					insertLengths.add(otArray[i]);
					hasInsertLength[otArray[i]] = true;
				}
			myFreq[0] += count;
			otArray = other.rgFlengths[otIdx];
			count = otArray[0];
			myFreq = this.rgFlFreq[myIdx];
			for (int i = 1; i <= count; i++)
				if (myFreq[otArray[i]]++ == 0 && !hasFragmentLength[otArray[i]]) {
					fragmentLengths.add(otArray[i]);
					hasFragmentLength[otArray[i]] = true;
				}
			myFreq[0] += count;
		}
	}

	private void mergeSampleArrays(FragmentLengthArrays other) {
		for (Map.Entry<String, Integer> e : smIndex.entrySet()) {
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
				if (myFreq[otArray[i]]++ == 0 && !hasInsertLength[otArray[i]]) {
					insertLengths.add(otArray[i]);
					hasInsertLength[otArray[i]] = true;
				}
			myFreq[0] += count;
			otArray = other.smFlengths[otIdx];
			count = otArray[0];
			myFreq = this.smFlFreq[myIdx];
			for (int i = 1; i <= count; i++)
				if (myFreq[otArray[i]]++ == 0 && !hasFragmentLength[otArray[i]]) {
					fragmentLengths.add(otArray[i]);
					hasFragmentLength[otArray[i]] = true;
				}
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
		this.size += other.size;
	}

	private void mergeFreqs(FragmentLengthFrequencies other, boolean arrayBased) {
		mergeSampleFreqs(other, arrayBased);
		mergeReadGroupFreqs(other, arrayBased);
	}

	private void mergeReadGroupFreqs(FragmentLengthFrequencies other,
			boolean arrayBased) {
		for (Map.Entry<String, Integer> e : rgIndex.entrySet()) {
			String rg = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.rgIndex.get(rg);
			if (otIdxI == null)
				continue;
			if (arrayBased) {
				arrayBasedMerge(other, myIdx, otIdxI, true, true);
				arrayBasedMerge(other, myIdx, otIdxI, true, false);
			} else {
				pgBasedMerge(other, myIdx, otIdxI, true, true);
				pgBasedMerge(other, myIdx, otIdxI, true, false);
			}
		}
	}

	private void mergeSampleFreqs(FragmentLengthFrequencies other,
			boolean arrayBased) {
		for (Map.Entry<String, Integer> e : smIndex.entrySet()) {
			String s = e.getKey();
			int myIdx = e.getValue();
			Integer otIdxI = other.smIndex.get(s);
			if (otIdxI == null)
				continue;
			if (arrayBased) {
				arrayBasedMerge(other, myIdx, otIdxI, false, true);
				arrayBasedMerge(other, myIdx, otIdxI, false, false);
			} else {
				pgBasedMerge(other, myIdx, otIdxI, false, true);
				pgBasedMerge(other, myIdx, otIdxI, false, false);
			}
		}
	}

	private void pgBasedMerge(FragmentLengthFrequencies ot, int myIdx,
			int otIdx, boolean rg, boolean il) {
		PriorityQueue<Integer> pq = il ? ot.insertLengths : ot.fragmentLengths;
		boolean[] hasLength = il ? hasInsertLength : hasFragmentLength;
		PriorityQueue<Integer> myPq = il ? insertLengths : fragmentLengths;

		Stack<Integer> st = new Stack<Integer>();
		st.ensureCapacity(pq.size());
		long[] myArray = rg ? (il ? rgIlFreq[myIdx] : rgFlFreq[myIdx])
				: (il ? smIlFreq[myIdx] : smFlFreq[myIdx]);
		long[] otArray = rg ? (il ? ot.rgIlFreq[otIdx] : ot.rgFlFreq[otIdx])
				: (il ? ot.smIlFreq[otIdx] : ot.smFlFreq[otIdx]);
		while (!pq.isEmpty()) {
			Integer next = pq.poll();
			myArray[next] += otArray[next];
			if (!hasLength[next]) {
				myPq.add(next);
				hasLength[next] = true;
			}
			st.push(next);
		}
		while (!st.isEmpty())
			pq.add(st.pop());
		myArray[0] += otArray[0];
	}

	private void arrayBasedMerge(FragmentLengthFrequencies ot, int myIdx,
			int otIdx, boolean rg, boolean il) {
		long[] myArray = rg ? (il ? rgIlFreq[myIdx] : rgFlFreq[myIdx])
				: (il ? smIlFreq[myIdx] : smFlFreq[myIdx]);
		long[] otArray = rg ? (il ? ot.rgIlFreq[otIdx] : ot.rgFlFreq[otIdx])
				: (il ? ot.smIlFreq[otIdx] : ot.smFlFreq[otIdx]);
		boolean[] hasLength = il ? hasInsertLength : hasFragmentLength;
		PriorityQueue<Integer> myPq = il ? insertLengths : fragmentLengths;
		for (int i = 0; i < otArray.length; i++) {
			myArray[i] += otArray[i];
			if (i > 0 && !hasLength[i]) {
				myPq.add(i);
				hasLength[i] = true;
			}
		}
	}

	private void increaseMaxLength(int newValue) {
		for (int i = 0; i < rgIlFreq.length; i++) {
			rgIlFreq[i] = Arrays.copyOf(rgIlFreq[i], newValue + 1);
			rgFlFreq[i] = Arrays.copyOf(rgFlFreq[i], newValue + 1);
		}
		for (int i = 0; i < smIlFreq.length; i++) {
			smIlFreq[i] = Arrays.copyOf(smIlFreq[i], newValue + 1);
			smIlFreq[i] = Arrays.copyOf(smIlFreq[i], newValue + 1);
		}
	}

	@Override
	public FragmentLengthSummary summary() {
		return new FragmentLengthSummary(this);
	}

}