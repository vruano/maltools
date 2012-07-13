package net.malariagen.gatk.utils;

import org.broadinstitute.sting.gatk.samples.Sample;

import net.sf.samtools.SAMReadGroupRecord;

public class ReadGroup {

	private SAMReadGroupRecord samRecord;
	
	private ReadGroupDB db;
	
	ReadGroup(ReadGroupDB db, SAMReadGroupRecord samRecord) {
		this.db = db;
		this.samRecord = samRecord;
	}
	
	public SAMReadGroupRecord getSAMRecord() {
		return samRecord;
	}
	
	public String getID() {
		return samRecord.getId();
	}
	
	public String getSampleID() {
		return samRecord.getSample();
	}
	
	public Sample getSample() {
		return db.sampleDb.getSample(samRecord);
	}
	
}
