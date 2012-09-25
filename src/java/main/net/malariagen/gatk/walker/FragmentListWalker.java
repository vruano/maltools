package net.malariagen.gatk.walker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import net.malariagen.utils.Nucleotide;
import net.malariagen.utils.io.TsvWriter;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

public class FragmentListWalker extends FragmentWalker<FragmentRecord,net.malariagen.gatk.walker.FragmentListWalker.Statistics> {

	@Output(shortName="o", fullName="output", doc="", required = true)
	protected File output;
	protected TsvWriter writer;
	
	@Argument(shortName="rg", fullName="outputReadGroup", doc="", required = false)
	protected boolean outputReadGroup = false;
	
	@Argument(shortName="sm", fullName="outputSample", doc="", required = false)
	protected boolean outputSample = false;
	
	@Argument(shortName="wf", fullName="outputFlags", doc="", required = false)
	protected boolean outputFlags = false;
	
	@Argument(shortName="mq", fullName="outputMappingQualities", doc="", required = false)
	protected boolean outputMappingQualities = false;
	
	@Argument(shortName="is", fullName="outputInsertStats", doc="", required= false)
	protected boolean outputInsertStats = false;

	@Argument(shortName="compress", fullName="compressOutput", doc="compress the output using GZIP. By default only output files finising with .gz are compressed", required = false)
	protected Boolean compressOutput = null;
	
	@Argument(shortName="readGC", fullName="outputReadGC", doc="GC content output is based on the reads as opposed to the fragment as imputed using the reference", required = false)
	protected boolean outputReadGC = false;
	
	@Argument(shortName = "mappedGC", fullName="outputMappedGC", doc="GC content output is based on mapped reads only. Implies -readGC when both reads are mapped.", required = false)
	protected boolean outputMappedGC = false;
	
	
	private int extrasCount = 0;
	private Object[] lineOutputs;
	
	
	@Override
	public void initialize() {
		super.initialize();
		if (compressOutput == null) 
			compressOutput = output.toString().matches(".*\\.gz") ? true : false;
		
		try {
			writer = compressOutput ? new TsvWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(output)))) : new TsvWriter(new FileWriter(output));
		} catch (IOException e) {
			throw new UserException("could not open output '" + output + "' to write",e);
		}
		try {
			writer.writeComment("FragmentListWalker v???");
			writer.writeComment(" maximumFragmentLength=" + this.maxLength);
			writer.writeComment(" minimumMappingQuality=" + this.minimumMappingQuality);
			writer.writeComment(" sort=" + this.order);
			writer.writeComment(" outputMappedGC=" + this.outputMappedGC);
			List<String> header = new ArrayList<String>(10);
			if (outputSample) header.add("SAMPLE");
			if (outputReadGroup) header.add("READ.GROUP");
			header.add("CHROM");
			header.add("START");
			header.add("END");
			if (outputFlags) header.add("FLAGS");
			if (outputMappingQualities) header.add("MQ.FIRST");
			if (outputMappingQualities) header.add("MQ.LAST");
			header.add("STRAND");
			header.add("LENGTH");
			if (outputInsertStats) {
			  header.add("I.START");
			  header.add("I.END");
			  header.add("I.LENGTH");
			}
			header.add("N.COUNT");
			header.add("GC.COUNT");
			if (outputReadGC) header.add("N.FIRST");
			if (outputReadGC) header.add("GC.FIRST");
			if (outputReadGC) header.add("N.LAST");
			if (outputReadGC) header.add("GC.LAST");
			writer.writeLine(header.toArray());
		} catch (IOException e) {
			throw new StingException("could not write in output file",e);
		}
		
		if (outputSample) extrasCount++;
		if (outputReadGroup) extrasCount++;
		if (outputFlags) extrasCount++;
		if (outputMappingQualities) extrasCount += 2;
		if (outputReadGC) extrasCount += 4;
		if (outputInsertStats) extrasCount += 3;
		
		lineOutputs = new Object[7 + extrasCount];
	}
	
	@Override
	public Statistics reduceFragmentInit() {
		return new Statistics();
	}

	@Override
	public FragmentRecord mapFragment(ReferenceContext ref,
			FragmentRecord fragment, FragmentMetaDataTracker metaDataTracker) {
		
		int nextIdx = 0;
		
		if (outputSample) lineOutputs[nextIdx++] = fragment.getFirstRead().getReadGroup().getSample();
		if (outputReadGroup) lineOutputs[nextIdx++] = fragment.getFirstRead().getReadGroup().getId();

		nextIdx = mapFragment$locInfo(fragment,lineOutputs,nextIdx);
		
		if (outputFlags) lineOutputs[nextIdx++] = FragmentFlag.flagsOf(fragment,this.minimumMappingQuality, this.maxLength);
		if (outputMappingQualities) nextIdx = mapFragment$mappingQualities(fragment,lineOutputs,nextIdx);
		lineOutputs[nextIdx++] = fragment.isProperFragment() ? (fragment.getPositiveStrand() ? "+" : "-") : "NA";
		if (fragment.isProperFragment()) {
		  lineOutputs[nextIdx++] = fragment.getLength();
		  if (outputInsertStats) {
		    lineOutputs[nextIdx++] = fragment.getInsertStart();
		    lineOutputs[nextIdx++] = fragment.getInsertEnd();
		    lineOutputs[nextIdx++] = fragment.getInsertLength();
		  }
		}
		else {
		  lineOutputs[nextIdx++] = "NA";
		  if (outputInsertStats) {
		    lineOutputs[nextIdx++] = "NA";
		    lineOutputs[nextIdx++] = "NA";
		    lineOutputs[nextIdx++] = "NA";
		  }
		}
		
		byte[] bases = outputMappedGC ? mapFragment$readBases(fragment) : (fragment.isProperFragment() ? ref.getBases() : new byte[0]);
		nextIdx = mapFragment$loadGCCounts(bases, lineOutputs, nextIdx);
		if (outputReadGC) {
		  nextIdx = mapFragment$loadGCCounts(fragment.getFirstReadReported().getReadBases(), lineOutputs, nextIdx);
		  nextIdx = mapFragment$loadGCCounts(fragment.getLastReadReported().getReadBases(), lineOutputs, nextIdx);
		}
		
		try {
			writer.writeLine(lineOutputs);
		} catch (IOException e) {
			throw new StingException("error writing in the output",e);
		}
		return fragment;
	}

	private int mapFragment$loadGCCounts(byte[] bases, Object[] lineOutputs, int nextIdx) {
		int gcCount = 0;
		int count = 0;
		
		for (byte b : bases) {
			switch(Nucleotide.fromByte(b)) {
			case N:
				break;
			case G:
			case C:
				gcCount++;
			default:
				count++;
			}
		}

		
		lineOutputs[nextIdx++] = count;
		lineOutputs[nextIdx++] = gcCount;
		return nextIdx;
	}
	

	
	private int mapFragment$locInfo(FragmentRecord fragment, Object[] dest, int idx) {
		SAMRecord refFirst = fragment.getReferenceFirstRead();
		SAMRecord refLast = fragment.getReferenceLastRead();
		if (refFirst.getReadUnmappedFlag() && refLast.getReadUnmappedFlag()) {
			dest[idx] = "NA"; dest[idx + 1] = "NA"; dest[idx + 2] = "NA";
		}
		else if (refFirst.getReadUnmappedFlag()) {
			dest[idx] = refLast.getReferenceName();
			dest[idx+1] = refLast.getAlignmentStart();
			dest[idx+2] = refLast.getAlignmentEnd();
		}
		else if (refLast.getReadUnmappedFlag()) {
			dest[idx] = refLast.getReferenceName();
			dest[idx + 1] = refFirst.getAlignmentStart();
			dest[idx + 2] = refFirst.getAlignmentEnd();
		}
		else {
			String firstName = refFirst.getReferenceName();
			String lastName = refLast.getReferenceName();
			if (firstName.equals(lastName))
				dest[idx] = firstName;
			else
				dest[idx] =  firstName + "/" + lastName;
			
			dest[idx + 1] = refFirst.getAlignmentStart();
			dest[idx + 2] = refLast.getAlignmentEnd();
		}
		return idx + 3;
	}
	
	private int mapFragment$mappingQualities(FragmentRecord fragment, Object[] dest, int idx) {
		SAMRecord refFirst = fragment.getReferenceFirstRead();
		SAMRecord refLast = fragment.getReferenceLastRead();
		dest[idx] = refFirst.getReadUnmappedFlag() ? "NA" : refFirst.getMappingQuality();
		dest[idx + 1] = refLast.getReadUnmappedFlag() ? "NA" : refLast.getMappingQuality();
		return idx + 2;
	}

	private byte[] mapFragment$readBases(FragmentRecord fragment) {
		SAMRecord first = fragment.getFirstRead();
		SAMRecord last = fragment.getLastRead();

		if (!outputMappedGC || (!first.getReadUnmappedFlag() && !last.getReadUnmappedFlag())) {
			  byte[] firstResult = fragment.getFirstRead().getReadBases();
			  byte[] secondResult = fragment.getLastRead().getReadBases();
			  byte[] result = Arrays.copyOf(firstResult, firstResult.length + secondResult.length);
			  System.arraycopy(secondResult, 0, result, firstResult.length, secondResult.length);
			  return result;			
		}
		
		if (first.getReadUnmappedFlag() && last.getReadUnmappedFlag()) {
		  return new byte[0];
		}
		else if (first.getReadUnmappedFlag()) {
		  return last.getReadBases();
		}
		else {
		  return first.getReadBases();
		}
	}


	public static class Statistics {
		
	}

	@Override
	public Statistics reduceFragment(FragmentRecord value, Statistics sum) {
		return sum;
	}
	
	@Override
	public void onFragmentTravesalDone(Statistics sum) {
		try {
			writer.close();
		} catch (IOException e) {
			throw new StingException("error closing output",e);
		}
	}

}
