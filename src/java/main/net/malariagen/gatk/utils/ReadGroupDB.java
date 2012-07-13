package net.malariagen.gatk.utils;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMReadGroupRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.samples.SampleDB;
import org.broadinstitute.sting.utils.exceptions.UserException;

public class ReadGroupDB {

	SampleDB sampleDb;
	Map<String,ReadGroup> readGroupsByID;
	Map<String,Set<ReadGroup>> readGroupsBySampleID;
	
	public ReadGroupDB(GenomeAnalysisEngine gae) {
		sampleDb = gae.getSampleDB();
		
		List<SAMReadGroupRecord> srgrList = gae.getReadsDataSource().getHeader().getReadGroups();
		readGroupsByID = new HashMap<String,ReadGroup>();
		readGroupsBySampleID = new HashMap<String,Set<ReadGroup>>();
		for (Sample s : sampleDb.getSamples()) { 
			if (readGroupsBySampleID.containsKey(s.getID()))
				throw new IllegalArgumentException("there is more than one sample with the same ID (" + s.getID() + ")");
			readGroupsBySampleID.put(s.getID(),new HashSet<ReadGroup>());
		}

		for (SAMReadGroupRecord srgr : srgrList) {
			ReadGroup rg = new ReadGroup(this,srgr);
			if (readGroupsByID.containsKey(srgr.getId()))
				throw new IllegalArgumentException("there is more than one read group with the same ID (" + srgr.getId() + ")");
			readGroupsByID.put(srgr.getId(),rg);
			String sampleID = rg.getSampleID();
			if (sampleID != null) {
				Set<ReadGroup> sampleRGs = readGroupsBySampleID.get(sampleID);
				if (sampleRGs == null)
					throw new IllegalArgumentException("read group " + rg.getID() + " belongs to an unknown sample " + sampleID);
				sampleRGs.add(rg);
			}
		}
		
		readGroupsByID = Collections.unmodifiableMap(readGroupsByID);
		for (Map.Entry<String,Set<ReadGroup>> srgs : readGroupsBySampleID.entrySet()) {
			Set<ReadGroup> rgs = srgs.getValue();
			if (rgs.size() == 0) {
				Logger.getLogger(getClass()).warn("there a sample (" + srgs.getKey() + ") without declared read-groups");
				rgs = Collections.emptySet();
			}
			if (rgs.size() != 1) {
				for (ReadGroup rg : rgs) {
					if (rg.getID().equals(srgs.getKey()))
						throw new UserException("a sample (" + rg.getID() + ") cannot have multiple read-groups if any of these has the same ID as the sample");
				}
				rgs = Collections.unmodifiableSet(rgs);
			}
			else {
				rgs = Collections.singleton(rgs.iterator().next());
			}
			srgs.setValue(rgs);
		}
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
	
	public int getReadGroupCount() {
		return readGroupsByID.size();
	}

	public ReadGroup getReadGroupByID(String gn) {
		return readGroupsByID.get(gn);
	}
	
	public Set<ReadGroup> getReadGroupsBySample(Sample s) {
		return readGroupsBySampleID.get(s.getID());
	}
	
	public Set<ReadGroup> getReadGroupsBySampleID(String sid) {
		return readGroupsBySampleID.get(sid);
	}
	
	public Set<ReadGroup> getReadGroupsBySampleOrReadGroupID(String sid) {
		if (readGroupsByID.containsKey(sid)) 
			return Collections.singleton(readGroupsByID.get(sid));
		else if (readGroupsBySampleID.containsKey(sid))
			return readGroupsBySampleID.get(sid);
		else
			return Collections.emptySet();
	}
}
