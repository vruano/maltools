package net.malariagen.gatk.math;

import java.util.LinkedHashMap;
import java.util.Map;

import net.malariagen.gatk.coverage.LocusCategory;


import com.google.gson.JsonObject;


public class IntegerSequenceDistributionSet {

	private Map<LocusCategory,ValueFrequencyArrayIntegerDistribution> categories;
	private ValueFrequencyArrayIntegerDistribution all;
	
	private IntegerSequenceDistributionSet() {
		
	}
	
	IntegerSequenceDistributionSet(IntegerSequenceCounterSet cs) {
		all = cs.all.toDistribution();
		categories = new LinkedHashMap<LocusCategory,ValueFrequencyArrayIntegerDistribution>(cs.categoriesByValue.size());
		for (Map.Entry<LocusCategory, IntegerCounter> e : cs.categoriesByValue.entrySet()) {
			categories.put(e.getKey(),e.getValue().toDistribution());
		}
	}
	
	public IntegerDistribution getCategoryDistribution(LocusCategory c) {
		return categories.get(c);
	}

	public IntegerDistribution getAllCategoriesDistribution() {
		return all;
	}
	
	
	public JsonObject toJsonObject() {
		JsonObject result = new JsonObject();
		result.add("ALL",all.toJsonObject());
		for (Map.Entry<LocusCategory, ValueFrequencyArrayIntegerDistribution> e : categories.entrySet()) 
			result.add(e.getKey().toString(),e.getValue().toJsonObject());
		return result;
	}
	
	public static IntegerSequenceDistributionSet fromJsonObject(JsonObject json) {
		IntegerSequenceDistributionSet result = new IntegerSequenceDistributionSet();
		result.loadJsonObject(json);
		return result;
	}

	private void loadJsonObject(JsonObject json) {
		all = ValueFrequencyArrayIntegerDistribution.fromJsonObject(json.get("ALL").getAsJsonObject());
		categories = new LinkedHashMap<LocusCategory,ValueFrequencyArrayIntegerDistribution>(LocusCategory.values().length);		
		
		for (LocusCategory c : LocusCategory.values())
			categories.put(c,ValueFrequencyArrayIntegerDistribution.fromJsonObject(json.get(c.toString()).getAsJsonObject()));
			
	}

	public static IntegerSequenceDistributionSet merge(IntegerSequenceDistributionSet ... sets) {
		if (sets.length == 0)
			throw new IllegalArgumentException("you need to specify at least one set");
		else if (sets.length == 1)
			return sets[0];
		IntegerSequenceDistributionSet result = new IntegerSequenceDistributionSet();
		ValueFrequencyArrayIntegerDistribution[] sequenceDistributions = new ValueFrequencyArrayIntegerDistribution[sets.length];
		for (int i = 0; i < sets.length; i++)
			sequenceDistributions[i] = sets[i].all;
		result.all = ValueFrequencyArrayIntegerDistribution.merge(sequenceDistributions);
		result.categories = new LinkedHashMap<LocusCategory,ValueFrequencyArrayIntegerDistribution>(sets[0].categories.size());
		for (Map.Entry<LocusCategory,ValueFrequencyArrayIntegerDistribution> e : sets[0].categories.entrySet()) {
			LocusCategory category = e.getKey();
			for (int i = 0; i < sets.length; i++)
				sequenceDistributions[i] = sets[i].categories.get(category);
			result.categories.put(category, ValueFrequencyArrayIntegerDistribution.merge(sequenceDistributions));
		}
		return result;
	}

}
