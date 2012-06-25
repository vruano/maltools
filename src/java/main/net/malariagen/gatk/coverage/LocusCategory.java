package net.malariagen.gatk.coverage;


import java.util.Collection;
import java.util.List;
import java.util.Map;

import net.malariagen.gatk.annotators.Constants;
import net.malariagen.gatk.gff.GFFFeature;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;

public enum LocusCategory {
	CODING, NON_CODING, GENIC, INTER_GENIC, FEATURELESS;

	public byte mask() {
		return (byte) (1 << this.ordinal());
	}

	public static byte CODING_MASK = CODING.mask();
	public static byte NON_CODING_MASK = NON_CODING.mask();
	public static byte GENIC_MASK = GENIC.mask();
	public static byte FEATURELESS_MASK = FEATURELESS.mask();
	public static byte INTER_GENIC_MASK = INTER_GENIC.mask();

	public static byte categoryMask(RefMetaDataTracker tracker) {
		List<RODRecordList> rrll = tracker.getBoundRodTracks();
		byte result = 0;
		for (RODRecordList rrl : rrll ) { 
			if (rrl.getName().equals(Constants.FEATURES_ROD_NAME)) 
				continue;
			for (GATKFeature ft : rrl) {
				Object o = ft.getUnderlyingObject();
				if (o instanceof GFFFeature) {
					switch (((GFFFeature) o).getType()) {
					case CDS:
						result |= CODING_MASK;
					case EXON:
					case GENE:
						result |= GENIC_MASK;
					}
				}
			}
		}
		if (result == 0)
			result |= FEATURELESS_MASK;
		if ((result & GENIC_MASK) == 0)
			result |= INTER_GENIC_MASK;
		if ((result & CODING_MASK) == 0)
			result |= NON_CODING_MASK;
		return result;
	}

	public static byte categoryMask(ReadMetaDataTracker tracker, SAMRecord read) {
		Map<Integer, Collection<GATKFeature>> tracks = tracker
				.getReadOffsetMapping(Constants.FEATURES_ROD_NAME);
		byte result = 0;
		int totalPositions = 0;
		int codingPositions = 0;
		int genicPositions = 0;
		int featurelessPositions = read.getReadLength() - tracks.size();
		for (Collection<GATKFeature> c : tracks.values()) {
			totalPositions++;
			if (c.size() == 0)
				featurelessPositions++;
			boolean isCoding = false;
			boolean isGenic = false;
			for (GATKFeature ft : c) {
				Object o = ft.getUnderlyingObject();
				if (o instanceof GFFFeature) {
					switch (((GFFFeature) o).getType()) {
					case CDS:
						isCoding = true;
					case EXON:
					case GENE:
						isGenic = true;
					}
				}
			}
			if (isCoding) codingPositions++;
			if (isGenic) genicPositions++;
		}
		int halfPositions = totalPositions >> 1;
		result |= (codingPositions >= halfPositions) ? CODING_MASK
				: NON_CODING_MASK;
		result |= (genicPositions >= halfPositions) ? GENIC_MASK
				: INTER_GENIC_MASK;
		if (featurelessPositions == totalPositions)
			result |= FEATURELESS_MASK;
		return result;
	}

}
