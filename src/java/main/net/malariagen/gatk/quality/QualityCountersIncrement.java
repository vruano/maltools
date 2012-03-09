package net.malariagen.gatk.quality;

import java.util.Arrays;

class QualityCountersIncrement {

	public QualityCountersIncrement() {
		this(200);
	}
	
	public QualityCountersIncrement(int maxReadLength) {
		baseQualities = new int[maxReadLength];
		baqs = new int[maxReadLength];
	}
	
	public void ensureCapacity(int l) {
		if (baseQualities.length < l) {
			baseQualities = Arrays.copyOf(baseQualities, l);
		    baqs = Arrays.copyOf(baqs, l);
		}
	}

	int sequence;
	
	int sample;
	
	int mappingQuality;
	
	int[] baseQualities;
	
	int baseQualityCount;
	
	int[] baqs;
	
	int baqCount;

	byte categories;

	public void applyTo(QualityCounters qs, boolean countSamples) {
		if (qs.mappingQuality != null) {
			qs.mappingQuality.addAllValue(mappingQuality,categories,sequence);
			qs.mappingQuality.addValue(mappingQuality, sample, categories, sequence);
		}
		if (qs.baq != null)
			for (int i = 0; i < baqCount; i++) {
				qs.baq.addAllValue(baqs[i],categories,sequence);
				qs.baq.addValue(baqs[i],sample,categories,sequence);
			}
		if (qs.baseQuality != null)
			for (int i = 0; i < baseQualityCount; i++) {
				qs.baseQuality.addAllValue(baseQualities[i],categories,sequence);
				qs.baseQuality.addValue(baseQualities[i],sample,categories,sequence);
			}
		
	}
	
}
