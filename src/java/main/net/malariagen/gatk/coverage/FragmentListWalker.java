package net.malariagen.gatk.coverage;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import net.malariagen.utils.Nucleotide;
import net.malariagen.utils.io.TsvWriter;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

public class FragmentListWalker extends FragmentWalker<FragmentRecord,net.malariagen.gatk.coverage.FragmentListWalker.Statistics> {

	@Output(shortName="o", fullName="output", doc="", required = true)
	protected File output;
	protected TsvWriter writer;
	
	@Argument(shortName="rg", fullName="outputReadGroup", doc="", required = false)
	protected boolean outputReadGroup = false;
	
	@Argument(shortName="sm", fullName="outputSample", doc="", required = false)
	protected boolean outputSample = false;

	@Argument(shortName="compress", fullName="compressOutput", doc="", required = false)
	protected Boolean compressOutput = null;
	
	
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
			List<String> header = new ArrayList<String>(10);
			if (outputSample) header.add("SAMPLE");
			if (outputReadGroup) header.add("READ.GROUP");
			header.add("CHROM");
			header.add("START");
			header.add("END");
			header.add("LENGTH");
			header.add("STRAND");
			header.add("I.START");
			header.add("I.END");
			header.add("I.LENGTH");
			header.add("N.COUNT");
			header.add("GC.COUNT");
			writer.writeLine(header.toArray());
		} catch (IOException e) {
			throw new StingException("could not write in output file",e);
		}
	}
	
	@Override
	public Statistics reduceFragmentInit() {
		return new Statistics();
	}

	@Override
	public FragmentRecord mapFragment(ReferenceContext ref,
			FragmentRecord fragment, FragmentMetaDataTracker metaDataTracker) {
		byte[] bases = ref.getBases();
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
		GenomeLoc loc = ref.getLocus();
		try {
			writer.writeLine(loc.getContig(),loc.getStart(),loc.getStop(),fragment.getLength(),
					fragment.getPositiveStrand() ? "+": "-",
							fragment.getInsertStart(),
							fragment.getInsertEnd(),fragment.getInsertLength(),count,gcCount);
		} catch (IOException e) {
			throw new StingException("error writing in the output",e);
		}
		return fragment;
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
