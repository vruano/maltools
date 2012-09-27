package net.malariagen.gatk.math;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;


public class IntegerCounterSet  {

	final Map<String,IntegerSampleCounterSet> samplesByName;
	private final IntegerSampleCounterSet[] samples;
	private final int sampleCount; 
	final IntegerSampleCounterSet all;
	
	public IntegerCounterSet(String[] sampleNames, String[] sequenceNames) {
		samplesByName = new HashMap<String,IntegerSampleCounterSet>(sampleNames.length);
		int nextIndex = 0;
		samples = new IntegerSampleCounterSet[sampleNames.length];
		for (String sampleName : sampleNames)	
			samplesByName.put(sampleName, samples[nextIndex++] = new IntegerSampleCounterSet(sequenceNames));
		sampleCount = samples.length;
		all = new IntegerSampleCounterSet(sequenceNames);
	}
	
	IntegerSampleCounterSet getSampleCounter(String name) {
		return samplesByName.get(name);
	}
	
	public void addAllValue(int value, byte categories, int sequence) {
		all.addValue(value, categories, sequence);
	}
	
	public void addSampleValues(int[] values, byte categories, int sequence, boolean useSum) {
		int total = 0;
		for (int i = 0; i < values.length; i++) {
			total += values[i];
			addValue(values[i],i,categories,sequence);
		}		
		if (useSum)
			addAllValue(total,categories,sequence);
	}
	
	public void addValue(int value, int sample, byte categories, int sequence) {
		samples[sample].addValue(value, categories, sequence);
	}

	public void applyCounterSet(IntegerCounterSet rhs) {
		all.applyCounterSet(rhs.all);
		for (int i = 0; i < sampleCount; i++) 
			samples[i].applyCounterSet(rhs.samples[i]);
	}

	public IntegerDistributionSet toCoverageDistributionSet() {
		return new IntegerDistributionSet(this);
	}
	
	public void printReport(PrintWriter pw) {
		pw.println("across all samples:");
		all.printReport("  ",pw);
		for (Map.Entry<String,IntegerSampleCounterSet> e : samplesByName.entrySet() ) {
			if (e.getValue().all.all.toDistribution().count() == 0) continue;			
			pw.println("Covarage in sample " + e.getKey() + ":");
			e.getValue().printReport("  ",pw);
		}
	}
}
