package net.malariagen.gatk.quality;

import java.io.IOException;
import java.io.Writer;

import net.malariagen.gatk.math.IntegerCounterSet;

class QualityCounters {

	public QualityCounters(boolean mq, boolean bq, boolean baq, String[] sampleNames, String[] sequenceNames) {
		mappingQuality = mq ? new IntegerCounterSet(sampleNames,sequenceNames) : null;
	    baseQuality = bq ? new IntegerCounterSet(sampleNames,sequenceNames) : null;
		this.baq = baq ? new IntegerCounterSet(sampleNames,sequenceNames) : null;
	}

	final IntegerCounterSet mappingQuality;

	final IntegerCounterSet baseQuality;

	final IntegerCounterSet baq;

	public void addFrom(QualityCounters other) {
		addFromCs(this.mappingQuality, other.mappingQuality);
		addFromCs(this.baseQuality, other.baseQuality);
		addFromCs(this.baq, other.baq);
	}

	private void addFromCs(IntegerCounterSet mine, IntegerCounterSet others) {
		if (mine != null && others != null)
			mine.applyCounterSet(others);
	}
	
	public void write(Writer m, Writer bq, Writer baq) throws IOException {
		write(mappingQuality,m);
		write(baseQuality,bq);
		write(this.baq,baq);
	}

	private void write(IntegerCounterSet ics, Writer w) throws IOException {
		if (ics == null || w == null)
			return;
		ics.toCoverageDistributionSet().write(w);
		
	}

}
