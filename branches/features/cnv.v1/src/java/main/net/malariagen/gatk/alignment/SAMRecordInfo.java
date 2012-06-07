package net.malariagen.gatk.alignment;

import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.concurrent.atomic.AtomicLong;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;

import net.malariagen.gatk.gff.GFFFeature;
import net.malariagen.gatk.uniqueness.UQNFeature;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

class SAMRecordInfo implements Comparable<SAMRecordInfo> {

	private final SAMRecord record;
	private Boolean isCoding = null;
	private String uniquenessString = null;
	private Integer uniquenessMedian = null;
	private String key;

	SAMRecordInfo(SAMRecord r) {
		record = r;
	}
	
	public String getKey() {
		return key;
	}

	public SAMRecordInfo(SAMRecord read, ReferenceContext ref,
			ReadMetaDataTracker metaDataTracker, boolean analyzeUniqueness,
			boolean analyzeLocusCategory) {
		this(read);
		if (analyzeUniqueness || analyzeLocusCategory) {
			List<GATKFeature> features = metaDataTracker.getAllCoveringRods();
			if (analyzeUniqueness)
				analyzeUniqueness(features);
			if (analyzeLocusCategory)
				analyzeLocusCategory(features);
		}
		key = read.getReadName() + "/" + (read.getFirstOfPairFlag() ? "1" : "2"); 
	}

	private void analyzeLocusCategory(List<GATKFeature> features) {
		int codingCount = 0;
		int totalCount = 0;
		for (GATKFeature f : features) {
			Object uo = f.getUnderlyingObject();
			if (!(uo instanceof GFFFeature))
				continue;
			totalCount++;
			GFFFeature gffFeature = (GFFFeature) uo;
			if (gffFeature.getType().isProteinCoding())
				codingCount++;
		}

		isCoding = codingCount > 0;
	}

	private void analyzeUniqueness(List<GATKFeature> features) {
		int[] uniquenessBytes = new int[record.getReadLength()];
		Arrays.fill(uniquenessBytes, 0);
		List<AlignmentBlock> ablocks = record.getAlignmentBlocks();
		ListIterator<AlignmentBlock> abIt = ablocks.listIterator();
		int scoreCount = 0;
		AlignmentBlock current = abIt.next();// (abIt.nextIndex() == ablocks.size()) ? abIt.previous() : abIt.next();
		int cbStart = current.getReferenceStart();
		int cbEnd = cbStart + current.getLength() - 1;
		for (GATKFeature f : features) {
			Object uo = f.getUnderlyingObject();
			if (!(uo instanceof UQNFeature))
				continue;
			UQNFeature uqn = (UQNFeature) uo;
			int score = uqn.getScore();
			int refPos = uqn.getStart();
			boolean foundIt = false;
			if (cbStart > refPos) { // forward search.
				while (cbStart > refPos && abIt.nextIndex() < ablocks.size()) {
					current = abIt.next();
					cbStart = current.getReferenceStart();
					cbEnd = cbStart + current.getLength()  - 1;
				}
				foundIt = cbStart <= refPos && cbEnd >= refPos;
			}
			else if (cbEnd < refPos) { // backward search.
				while (cbEnd < refPos && abIt.previousIndex() >= 0) {
					current = abIt.previous();
					cbStart = current.getReferenceStart();
					cbEnd = cbStart + current.getLength() - 1;
				}
				foundIt = cbStart <= refPos && cbEnd >= refPos;
			}
			else {
				foundIt = true;
			}
			if (foundIt) {
				int offset = current.getReadStart() + (refPos - cbStart) - 1;
				uniquenessBytes[offset] = score;
			    scoreCount++;
			}
		}
		if (scoreCount == 0) { 
			uniquenessString = null;
		    uniquenessMedian = null;
		    return;
		}
		uniquenessString = uniquenessArrayToString(uniquenessBytes);
		Arrays.sort(uniquenessBytes);
		int start = Arrays.binarySearch(uniquenessBytes, 1);
		if (start > -uniquenessBytes.length - 1) {
			if (start < 0)
				start = -start - 1;
			int size = uniquenessBytes.length -start;
			if (size == 0) 
				uniquenessMedian = null;
			else if (size == 1)
				uniquenessMedian = uniquenessBytes[start];
			else if (size == 2)
				uniquenessMedian = (uniquenessBytes[start] + uniquenessBytes[start + 1]) >> 1;
			else {
			  int midSize = size >> 1;
			  boolean even = ((size & 1) == 0);
			  int floorIndex = start + midSize + (even ? 0 : 1);
			  int ceilIndex = floorIndex + (even  ? 1 : 0);
			  if (ceilIndex >= uniquenessBytes.length) 
				throw new IllegalStateException(" " + start + " " + even + " " + size + " " + uniquenessBytes[start] + " " + uniquenessBytes[start - 1]);
			  uniquenessMedian = ((uniquenessBytes[floorIndex] + uniquenessBytes[ceilIndex]) >> 1);
			}
	
		} else
			uniquenessMedian = null;
	}

	private String uniquenessArrayToString(int[] uniquenessScores) {
		StringBuffer bf = new StringBuffer(uniquenessScores.length);
		int i = 0;
		while (i < uniquenessScores.length) {
			int score = uniquenessScores[i];
			int first = i;
			while (i < uniquenessScores.length && uniquenessScores[i] == score)
				i++;
			if (score == 0) {
				bf.append(- (i - first));
			}
			else {
				bf.append( (i - first )).append('*').append(score);
			}
			if (i < uniquenessScores.length) bf.append('+');
		}
		return bf.toString();
	}

	public SAMRecord getRecord() {
		return record;
	}


	public String getUniquenessScoreString() {
		return uniquenessString;
	}

	public String getUniquenessMedian() {
		return uniquenessMedian == null ? null : uniquenessMedian.toString();
	}

	public String getLocusCategory() {
		return isCoding == null ? null : ((isCoding) ? "1" : "0");
	}

	String getField(SAMRecordInfoField f) {
		String v = null;
		switch (f) {
		case SQ:
			v = record.getReadString();
			break;
		case BQ:
			v = record.getBaseQualityString();
			break;
		case MQ:
			v = "" + record.getMappingQuality();
			break;
		case CG:
			v = "" + record.getCigarString();
			break;
		case UQ:
			v = getUniquenessScoreString();
			break;
		case UM:
			v = getUniquenessMedian();
			break;
		case SP:
			v = "" + record.getAlignmentStart();
			break;
		case EP:
			v = "" + record.getAlignmentEnd();
			break;
		case LC:
			v = getLocusCategory();
			break;
		//case FL:
		//	v = "" + record.getFlags();
		}
		return v;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o)
			return true;
		if (o == null)
			return false;
		if (!(o instanceof SAMRecordInfo))
			return false;
		return equals((SAMRecordInfo)o);
	}
	
	public boolean equals(SAMRecordInfo other) {
		if (this == other)
			return true;
		if (other == null)
			return false;
		
		if (!this.getKey().equals(other.getKey()))
			return false;
		
		for (SAMRecordInfoField f : SAMRecordInfoField.values()) {
			String thisValue = this.getField(f);
			String otherValue = other.getField(f);
			if (thisValue == otherValue) continue;
			if (thisValue == null && otherValue != null) return false;
			if (thisValue != null && otherValue == null) return false;
			if (!thisValue.equals(otherValue)) return false;			
		}
		return true;
	}

	public boolean equalsAnalysis(SAMRecordInfo other) {
		if (this == other)
			return true;
		if (other == null)
			return false;
		
		if (!this.getKey().equals(other.getKey())) {
			return false;
		}
		
		for (SAMRecordInfoField f : SAMRecordInfoField.values()) {
			String thisValue = this.getField(f);
			String otherValue = other.getField(f);
			if (!thisValue.equals(otherValue)) return false;
		}
		return true;
	}
	
	
	@Override
	public int compareTo(SAMRecordInfo r2) {
		
		SAMRecordInfo r1 = this;
		if (r1.equals(r2))
			return 0;
		int contig1 = r1.getRecord().getReferenceIndex();
		int contig2 = r2.getRecord().getReferenceIndex();
		if (contig1 < contig2)
			return -1;
		if (contig2 < contig1)
			return 1;
		int start1 = r1.getRecord().getAlignmentStart();
		int start2 = r2.getRecord().getAlignmentStart();
		if (start1 < start2)
			return -1;
		else if (start1 > start2)
			return 1;
		else {
			int end1 = r1.getRecord().getAlignmentEnd();
			int end2 = r2.getRecord().getAlignmentEnd();
			if (end1 < end2)
				return -1;
			else if (end1 > end2)
				return 1;
			else {
				long serial1 = r1.getSerialNumber();
				long serial2 = r2.getSerialNumber();
				if (serial1 < serial2)
					return -1;
				else if (serial2 < serial1)
					return 1;
				else {
					throw new RuntimeException("should never happend " + serial1 + " == " + serial2);
				}
			}
		}
	}

	private final static AtomicLong nextSerialNumber = new AtomicLong(Long.MIN_VALUE);

	private long serialNumber = nextSerialNumber.getAndIncrement();

	public long getSerialNumber() {
		return serialNumber;
	}

}
