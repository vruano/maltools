package net.malariagen.gatk.math;

import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.LinkedHashMap;
import java.util.Map;

import net.malariagen.gatk.coverage.LocusCategory;


import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;


public class IntegerDistributionSet {

	static {
		forceClassLoading();
	}
	
	private Map<String,IntegerSampleDistributionSet> samples;
	private IntegerSampleDistributionSet all;
	
	public IntegerDistributionSet() {
		
	}
	
	// this is a fix for a very weird bug when running the stuff in Sanger's farm.. that results in a Zip
	private static void forceClassLoading() {
		// ValueFrequencyArrayIntegerDistribution v = new ValueFrequencyArrayIntegerDistribution(20,10 + 20,10 + 4 * 10,new int[] { 1,2 }, new int[] { 10,10 });
		
	}

	IntegerDistributionSet(IntegerCounterSet cs) {
		all = new IntegerSampleDistributionSet(cs.all);
		samples = new LinkedHashMap<String,IntegerSampleDistributionSet>(cs.samplesByName.size());
		for (Map.Entry<String, IntegerSampleCounterSet> e : cs.samplesByName.entrySet()) {
			samples.put(e.getKey(),new IntegerSampleDistributionSet(e.getValue()));
		}
	}
	
	public IntegerSampleDistributionSet getSampleDistributionSet(String seq) {
		return samples.get(seq);
	}

	public IntegerSampleDistributionSet getAllSamplesDistributionSet() {
		return all;
	}
	
	
	public JsonObject toJsonObject() {
		JsonObject result = new JsonObject();
		result.add("ALL",all.toJsonObject());
		for (Map.Entry<String,IntegerSampleDistributionSet> e : samples.entrySet()) 
			result.add(e.getKey(),e.getValue().toJsonObject());
		return result;
	}
	
	public static IntegerDistributionSet fromJsonObject(JsonObject json) {
		IntegerDistributionSet result = new IntegerDistributionSet();
		result.loadJsonObject(json);
		return result;
	}

	private void loadJsonObject(JsonObject json) {
		all = IntegerSampleDistributionSet.fromJsonObject(json.get("ALL").getAsJsonObject());
		samples = new LinkedHashMap<String,IntegerSampleDistributionSet>(LocusCategory.values().length);		
	    for ( Map.Entry<String,JsonElement> o : json.entrySet() ) {
	    	String sampleName = o.getKey();
	    	if (sampleName.equalsIgnoreCase("ALL")) 
	    		all = IntegerSampleDistributionSet.fromJsonObject(o.getValue().getAsJsonObject());
	    	else 
	    		samples.put(sampleName,IntegerSampleDistributionSet.fromJsonObject(o.getValue().getAsJsonObject()));
	    }
	}
	
	public static IntegerDistributionSet merge(IntegerDistributionSet ... sets) {
		if (sets.length == 0)
			throw new IllegalArgumentException("you need to specify at least one set");
		else if (sets.length == 1)
			return sets[0];
		IntegerDistributionSet result = new IntegerDistributionSet();
		IntegerSampleDistributionSet[] sampleDistributions = new IntegerSampleDistributionSet[sets.length];
		for (int i = 0; i < sets.length; i++)
			sampleDistributions[i] = sets[i].all;
		result.all = IntegerSampleDistributionSet.merge(sampleDistributions);
		result.samples = new LinkedHashMap<String,IntegerSampleDistributionSet>(sets[0].samples.size());
		for (Map.Entry<String,IntegerSampleDistributionSet> e : sets[0].samples.entrySet()) {
			String sequence = e.getKey();
			for (int i = 0; i < sets.length; i++)
				sampleDistributions[i] = sets[i].samples.get(sequence);
			result.samples.put(sequence, IntegerSampleDistributionSet.merge(sampleDistributions));
		}
		return result;
	}
	
	public void write(Writer w) throws IOException {
		w.write(toJsonObject().toString());
		w.flush();
	}
	
	public static IntegerDistributionSet read(Reader r) throws IOException {
		return fromJsonObject(new JsonParser().parse(r).getAsJsonObject());
	}

}
