package net.malariagen.gatk.math;

import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.LinkedHashMap;
import java.util.Map;


import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

public class IntegerSampleDistributionSet {

	private Map<String, IntegerSequenceDistributionSet> sequences;
	private IntegerSequenceDistributionSet all;

	private static JsonParser JSON_PARSER = new JsonParser();

	private IntegerSampleDistributionSet() {

	}

	IntegerSampleDistributionSet(IntegerSampleCounterSet cs) {
		all = new IntegerSequenceDistributionSet(cs.all);
		sequences = new LinkedHashMap<String, IntegerSequenceDistributionSet>(
				cs.sequencesByName.size());
		for (Map.Entry<String, IntegerSequenceCounterSet> e : cs.sequencesByName
				.entrySet()) {
			sequences.put(e.getKey(),
					new IntegerSequenceDistributionSet(e.getValue()));
		}
	}

	public IntegerSequenceDistributionSet getSequenceDistributionSet(String seq) {
		return sequences.get(seq);
	}

	public IntegerSequenceDistributionSet getAllSequencesDistributionSet() {
		return all;
	}

	public JsonObject toJsonObject() {
		JsonObject result = new JsonObject();
		result.add("ALL", all.toJsonObject());
		for (Map.Entry<String, IntegerSequenceDistributionSet> e : sequences
				.entrySet())
			result.add(e.getKey(), e.getValue().toJsonObject());
		return result;
	}

	public static IntegerSampleDistributionSet read(Reader r)
			throws IOException {

		IntegerSampleDistributionSet result = fromJsonObject(JSON_PARSER.parse(r).getAsJsonObject());
		if (result == null) {
			throw new IOException();
		}
		return result;
	}

	public void write(Writer w) throws IOException {
		w.write(this.toJsonObject().toString());
	}

	public static IntegerSampleDistributionSet fromJsonObject(JsonObject json) {
		IntegerSampleDistributionSet result = new IntegerSampleDistributionSet();
		result.loadJsonObject(json);
		return result;
	}

	private void loadJsonObject(JsonObject json) {
		all = IntegerSequenceDistributionSet.fromJsonObject(json
				.get("ALL").getAsJsonObject());
		sequences = new LinkedHashMap<String, IntegerSequenceDistributionSet>(20);
	    for ( Map.Entry<String,JsonElement> o : json.entrySet() ) {
	    	String sequenceName = o.getKey();
	    	if (sequenceName.equalsIgnoreCase("ALL")) 
	    		all = IntegerSequenceDistributionSet.fromJsonObject(o.getValue().getAsJsonObject());
	    	else 
	    		sequences.put(sequenceName,IntegerSequenceDistributionSet.fromJsonObject(o.getValue().getAsJsonObject()));
	    }
	}
	
	public static IntegerSampleDistributionSet merge(IntegerSampleDistributionSet ... sets) {
		if (sets.length == 0)
			throw new IllegalArgumentException("you need to specify at least one set");
		else if (sets.length == 1)
			return sets[0];
		IntegerSampleDistributionSet result = new IntegerSampleDistributionSet();
		IntegerSequenceDistributionSet[] sequenceDistributions = new IntegerSequenceDistributionSet[sets.length];
		for (int i = 0; i < sets.length; i++)
			sequenceDistributions[i] = sets[i].all;
		result.all = IntegerSequenceDistributionSet.merge(sequenceDistributions);
		result.sequences = new LinkedHashMap<String,IntegerSequenceDistributionSet>(sets[0].sequences.size());
		for (Map.Entry<String,IntegerSequenceDistributionSet> e : sets[0].sequences.entrySet()) {
			String sequence = e.getKey();
			for (int i = 0; i < sets.length; i++)
				sequenceDistributions[i] = sets[i].sequences.get(sequence);
			result.sequences.put(sequence, IntegerSequenceDistributionSet.merge(sequenceDistributions));
		}
		return result;
	}

}
