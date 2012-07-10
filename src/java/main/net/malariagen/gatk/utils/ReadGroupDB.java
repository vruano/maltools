package net.malariagen.gatk.utils;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMReadGroupRecord;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.samples.SampleDB;

public class ReadGroupDB {

	SampleDB sampleDb;
	Map<String,ReadGroup> readGroupsByID;
	public ReadGroupDB(GenomeAnalysisEngine gae) {
		sampleDb = gae.getSampleDB();
		
		List<SAMReadGroupRecord> srgrList = gae.getReadsDataSource().getHeader().getReadGroups();
		readGroupsByID = new HashMap<String,ReadGroup>();
		for (SAMReadGroupRecord srgr : srgrList)
			readGroupsByID.put(srgr.getId(),new ReadGroup(this,srgr));
		
		readGroupsByID = Collections.unmodifiableMap(readGroupsByID);
	}
	
	public Map<String,ReadGroup> getReadGroupsByID() {
		return readGroupsByID;
	}
	
	public Set<String> getReadGroupIDs() {
		return readGroupsByID.keySet();
	}
	
	public SampleDB getSampleDB() {
		return sampleDb;
	}
	
	public int getGroupCount() {
		return readGroupsByID.size();
	}

	public ReadGroup getReadGroupByID(String gn) {
		return readGroupsByID.get(gn);
	}
	
	
	
}
