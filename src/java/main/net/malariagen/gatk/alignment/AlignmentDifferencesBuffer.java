package net.malariagen.gatk.alignment;

import java.io.PrintWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.utils.collections.Pair;

import net.sf.samtools.SAMRecord;

public class AlignmentDifferencesBuffer {

	private int bufferSize;
	private int windowSize;

	private StringBuffer buffer = new StringBuffer();

	private SortedSet<Pair<SAMRecordInfo, SAMRecordInfo>> matchedSet;
	private SortedSet<SAMRecordInfo> unmatchedSet;
	private Map<String, SAMRecordInfo> unmatchedByName;

	private SAMReaderID leftReaderID;

	private PrintWriter output;
	private int contigIndex = -1;
	private static SAMRecordInfoField[] INFO_FIELDS;
	private boolean[] differs;
	private String[] leftValues;
	private String[] rightValues;
	private CompareAlignmentsWalker walker;
	private int usedSize;
	private static final Comparator<? super Pair<SAMRecordInfo, SAMRecordInfo>> PAIR_COMPARATOR_BY_POSITION = new Comparator<Pair<SAMRecordInfo, SAMRecordInfo>>() {
		public int compare(Pair<SAMRecordInfo, SAMRecordInfo> p1,
				Pair<SAMRecordInfo, SAMRecordInfo> p2) {

			int start11 = p1.first.getRecord().getAlignmentStart();
			int start12 = p1.second.getRecord().getAlignmentStart();
			int start1 = start11 < start12 ? start11 : start12;
			int start21 = p2.first.getRecord().getAlignmentStart();
			int start22 = p2.second.getRecord().getAlignmentStart();
			int start2 = start21 < start22 ? start21 : start22;
			if (start1 < start2)
				return -1;
			if (start2 < start1)
				return 1;
			int end11 = p1.first.getRecord().getAlignmentEnd();
			int end12 = p1.second.getRecord().getAlignmentEnd();
			int end1 = end11 < end12 ? end12 : end11;
			int end21 = p2.first.getRecord().getAlignmentEnd();
			int end22 = p2.second.getRecord().getAlignmentEnd();
			int end2 = end21 < end22 ? end22 : end21;
			if (end1 < end2)
				return -1;
			if (end1 > end2)
				return 1;
			else
				return p1.first.compareTo(p2.first);
		}
	};

	public AlignmentDifferencesBuffer(CompareAlignmentsWalker walker) {
		bufferSize = walker.bufferSize;
		if (bufferSize < 0)
			bufferSize = Integer.MAX_VALUE;
		usedSize = 0;
		windowSize = walker.windowSize;
		if (windowSize < 0)
			windowSize = Integer.MAX_VALUE;
		output = walker.output;
		this.walker = walker;
		leftReaderID = walker.leftReaderID;
		unmatchedSet = new TreeSet<SAMRecordInfo>();
		matchedSet = new TreeSet<Pair<SAMRecordInfo, SAMRecordInfo>>(
				PAIR_COMPARATOR_BY_POSITION);
		unmatchedByName = new HashMap<String, SAMRecordInfo>(bufferSize);
	}

	public void add(SAMRecordInfo r) {
		int rContigIndex = r.getRecord().getReferenceIndex();
		if (contigIndex >= 0 && contigIndex != rContigIndex)
			flush();
		contigIndex = rContigIndex;
		SAMRecordInfo s = unmatchedByName.get(r.getKey());
		if (s != null) {
			unmatchedSet.remove(s);
			unmatchedByName.remove(s.getKey());
			int rsComp = r.compareTo(s);
			if (rsComp == 0) {
				usedSize--;
				return; // sure that buffer max size won't be exceeded.
			}
			if (r.compareTo(s) < 0)
				matchedSet.add(new Pair<SAMRecordInfo, SAMRecordInfo>(r, s));
			else
				matchedSet.add(new Pair<SAMRecordInfo, SAMRecordInfo>(s, r));
			usedSize++;
		} else {
			unmatchedSet.add(r);
			unmatchedByName.put(r.getKey(), r);
			usedSize++;
		}
		SAMRecordInfo firstUnmatched = unmatchedSet.isEmpty() ? null
				: unmatchedSet.first();
		Pair<SAMRecordInfo, SAMRecordInfo> firstMatched = matchedSet.isEmpty() ? null
				: matchedSet.first();
		int matchedStart = firstMatched == null ? Integer.MAX_VALUE
				: firstMatched.first.getRecord().getAlignmentStart();
		int unmatchedStart = firstUnmatched == null ? Integer.MAX_VALUE
				: firstUnmatched.getRecord().getAlignmentStart();
		int lastPosition = r.getRecord().getAlignmentStart();
		int firstPosition = matchedStart < unmatchedStart ? matchedStart
				: unmatchedStart;
		while (usedSize > bufferSize
				|| (lastPosition - firstPosition >= windowSize)) {
			if (matchedStart < unmatchedStart) {
				processPair(firstMatched.first, firstMatched.second);
				matchedSet.remove(firstMatched);
				usedSize -= 2;
				firstMatched = matchedSet.isEmpty() ? null : matchedSet.first();
				matchedStart = firstMatched == null ? Integer.MAX_VALUE
						: firstMatched.first.getRecord().getAlignmentStart();
			} else {
				processSingleton(firstUnmatched);
				unmatchedSet.remove(firstUnmatched);
				unmatchedByName.remove(firstUnmatched.getKey());
				firstUnmatched = unmatchedSet.isEmpty() ? null : unmatchedSet
						.first();
				unmatchedStart = firstUnmatched == null ? Integer.MAX_VALUE
						: firstUnmatched.getRecord().getAlignmentStart();
				usedSize--;
			}
			firstPosition = matchedStart < unmatchedStart ? matchedStart
					: unmatchedStart;
		}
	}

	public void flush() {
		SAMRecordInfo firstUnmatched = unmatchedSet.isEmpty() ? null
				: unmatchedSet.first();
		Pair<SAMRecordInfo, SAMRecordInfo> firstMatched = matchedSet.isEmpty() ? null
				: matchedSet.first();
		int matchedStart = firstMatched == null ? Integer.MAX_VALUE
				: firstMatched.first.getRecord().getAlignmentStart();
		int unmatchedStart = firstUnmatched == null ? Integer.MAX_VALUE
				: firstUnmatched.getRecord().getAlignmentStart();
		while (usedSize > 0) {
			if (matchedStart < unmatchedStart) {
				processPair(firstMatched.first, firstMatched.second);
				matchedSet.remove(firstMatched);
				usedSize -= 2;
				firstMatched = matchedSet.isEmpty() ? null : matchedSet.first();
				matchedStart = firstMatched == null ? Integer.MAX_VALUE
						: firstMatched.first.getRecord().getAlignmentStart();
			} else {
				processSingleton(firstUnmatched);
				unmatchedSet.remove(firstUnmatched);
				unmatchedByName
						.remove(firstUnmatched.getKey());
				firstUnmatched = unmatchedSet.isEmpty() ? null : unmatchedSet
						.first();
				unmatchedStart = firstUnmatched == null ? Integer.MAX_VALUE
						: firstUnmatched.getRecord().getAlignmentStart();
				usedSize--;
			}
		}
		// some checkings.
		if (!unmatchedSet.isEmpty() || !matchedSet.isEmpty())
			throw new RuntimeException("bug!!!");
		contigIndex = -1;
	}

	private void processSingleton(SAMRecordInfo r) {
		buffer.setLength(0);
		SAMRecord read = r.getRecord();
		buffer.append(read.getReferenceName()).append('\t');
		buffer.append(read.getAlignmentStart()).append('\t');
		buffer.append(read.getAlignmentEnd()).append('\t');
		buffer.append(read.getReadName());
		if (walker.getReaderIDFromRecord(read).equals(leftReaderID))
			buffer.append('\t').append('<');
		else
			buffer.append('\t').append('>');

		for (SAMRecordInfoField f : SAMRecordInfoField.values()) {
			String v = r.getField(f);
			if (f == SAMRecordInfoField.SP || f == SAMRecordInfoField.EP)
				continue; // not start and end position.
			if (v != null)
				buffer.append('\t').append(f).append(':').append(v);
		}
		output.println(buffer);
	}

	private void processPair(SAMRecordInfo r, SAMRecordInfo s) {
		if (r.equals(s))
			return;
		SAMReaderID rReaderID = walker.getReaderIDFromRecord(r.getRecord());
		SAMRecordInfo left, right;
		if (rReaderID.equals(leftReaderID)) {
			left = r;
			right = s;
		} else {
			left = s;
			right = r;
		}
		outputPair(left, right);
	}

	private void outputPair(SAMRecordInfo left, SAMRecordInfo right) {
		buffer.setLength(0);
		SAMRecord read1 = left.getRecord();
		SAMRecord read2 = left.getRecord();
		int start1 = read1.getAlignmentStart();
		int start2 = read2.getAlignmentStart();
		int end1 = read1.getAlignmentEnd();
		int end2 = read2.getAlignmentEnd();
		int start = start1 < start2 ? start1 : start2;
		int end = end1 < end2 ? end2 : end1;
		buffer.append(read1.getReferenceName()).append('\t');
		buffer.append(start).append('\t');
		buffer.append(end).append('\t');
		buffer.append(read1.getReadName());

		INFO_FIELDS = SAMRecordInfoField.values();
		differs = new boolean[INFO_FIELDS.length];
		int diffCount = 0;
		int leftDiffCount = 0;
		int rightDiffCount = 0;
		int sameCount = 0;
		leftValues = new String[INFO_FIELDS.length];
		rightValues = new String[INFO_FIELDS.length];
		for (int i = 0; i < INFO_FIELDS.length; i++) {
			leftValues[i] = left.getField(INFO_FIELDS[i]);
			rightValues[i] = right.getField(INFO_FIELDS[i]);
			differs[i] = !(leftValues[i] != null ? leftValues[i]
					.equals(rightValues[i]) : rightValues[i] == null);
			if (differs[i]) {
				diffCount++;
				if (leftValues[i] != null)
					leftDiffCount++;
				if (rightValues[i] != null)
					rightDiffCount++;
			}
			if (!differs[i] && leftValues[i] != null)
				sameCount++;
		}
		
		// first we output common stuff behind '='.
		if (sameCount > 0) {
			buffer.append('\t').append('=');
			for (int i = 0; i < INFO_FIELDS.length; i++) {
				if (differs[i])
					continue;
				if (INFO_FIELDS[i] == SAMRecordInfoField.SP
						|| INFO_FIELDS[i] == SAMRecordInfoField.EP)
					continue; // not report common values for start and end
								// position.
				if (leftValues[i] == null)
					continue;
				buffer.append('\t').append(INFO_FIELDS[i]).append(':')
						.append(leftValues[i]);
			}
		}
		if (diffCount > 0) {
			if (leftDiffCount > 0) {
				buffer.append('\t').append('<');
				for (int i = 0; i < INFO_FIELDS.length; i++) {
					if (!differs[i])
						continue;
					if (leftValues[i] == null)
						continue;
					buffer.append('\t').append(INFO_FIELDS[i]).append(':')
							.append(leftValues[i]);
				}
			}
			if (rightDiffCount > 0) {
				buffer.append('\t').append('>');
				for (int i = 0; i < INFO_FIELDS.length; i++) {
					if (!differs[i])
						continue;
					if (rightValues[i] == null)
						continue;
					buffer.append('\t').append(INFO_FIELDS[i]).append(':')
							.append(rightValues[i]);
				}
			}
		}
		output.println(buffer);
	}


}
