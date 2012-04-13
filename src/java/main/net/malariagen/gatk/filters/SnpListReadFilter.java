package net.malariagen.gatk.filters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.malariagen.gatk.csl.CSLFeature;
import net.malariagen.utils.NucleotideIUPAC;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

public class SnpListReadFilter extends ReadFilter {

	private boolean cslProvided;
	
	private byte[][] refMask;
	
	public SnpListReadFilter() {
	}
	
	public SnpListReadFilter(GenomeAnalysisEngine engine) {
		this();
		initialize(engine);
	}
	
	@Override
	public void initialize(GenomeAnalysisEngine engine) {
		List<ReferenceOrderedDataSource> sources = engine.getRodDataSources();
		ReferenceDataSource refSource = engine.getReferenceDataSource();
		IndexedFastaSequenceFile reference = refSource.getReference();
		GenomeLocParser glp = new GenomeLocParser(reference);
		
		List<SAMSequenceRecord> seqRecords = reference.getSequenceDictionary().getSequences();
		refMask = new byte[seqRecords.size()][];
		int nextIdx = 0;
		for (SAMSequenceRecord seqRecord : seqRecords) 
			refMask[nextIdx++] = reference.getSequence(seqRecord.getSequenceName()).getBases();
		Logger logger = Logger.getLogger(SnpListReadFilter.class);
		for (ReferenceOrderedDataSource source : sources) {
			if (source.getName().equals("csl")) {
				if (!cslProvided) {
					cslProvided = true;
					logger.info("CSL: creating reference mask from snp-list file ...");
				}
				parseCslRod(seqRecords,source,glp);
			}
		}
		if (cslProvided)
			logger.info("CSL: finished creating the mask");
	}
	
	private void parseCslRod(List<SAMSequenceRecord> seqs, ReferenceOrderedDataSource source, GenomeLocParser glp) {
		for (int i = 0; i < seqs.size(); i++) {
			SAMSequenceRecord seq = seqs.get(i);
			byte[] mask = refMask[i];
			GenomeLoc gl = glp.createGenomeLoc(seq.getSequenceName(), 1, seq.getSequenceLength());
			LocationAwareSeekableRODIterator it = source.seek(gl);
			while (it.hasNext()) { 
				RODRecordList rl = it.next();
				for (GATKFeature feature : rl) {
					Object ulo = feature.getUnderlyingObject();
					if (ulo instanceof CSLFeature) {
						CSLFeature csl = (CSLFeature) ulo;
						mask[csl.getStart() - 1] = csl.getAlternative();
					}
					else {
						int start = feature.getStart();
						int end = feature.getEnd();
						if (start == end)
							mask[start - 1] = 'N';
						else
							for (int j = start; j <= end; j++)
								mask[j - 1] = 'N';
					}

				}	
			}
				
		}
	}

	
	public boolean filterOut(GenomeLoc gl) {
		if (!cslProvided)
			return false;
		
		int refIndex = gl.getContigIndex();
	    byte[] mask = refMask[refIndex];
	    int start = gl.getStart();
	    int stop = gl.getStop();
	    if (start == stop)
	      return mask[start - 1] == 0;
	    else {
	      for (int i = start; i <= stop; i++) 
	    	  if (mask[i - 1] != 0) return true;
	      return false;
	    }
	}
	
	@Override
	public boolean filterOut(SAMRecord read) {
		if (!cslProvided)
			return false;
		if (read.getReadUnmappedFlag())
			return false;
		int refIndex = read.getReferenceIndex();
		byte[] mask = refMask[refIndex];
		
		byte[] bases = read.getReadBases();
		for (AlignmentBlock ab : read.getAlignmentBlocks()) {
			int refStart = ab.getReferenceStart() - 1;
			int abLen = ab.getLength();
			int readIdx = ab.getReadStart();
			for (int baseIdx = refStart; abLen > 0; abLen--) {
				if (NucleotideIUPAC.areCompatible(mask[baseIdx++],bases[readIdx++])) 
					return true;
			}
			
		}
		return false;
	}
}
