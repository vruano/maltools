package net.malariagen.gatk.math;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import net.malariagen.gatk.coverage.LocusCategory;


public class IntegerSequenceCounterSet {

	final Map<LocusCategory,IntegerCounter> categoriesByValue;
	private final IntegerCounter[] categories;
	final IntegerCounter all;
	private final int categoryCount;
	
	IntegerSequenceCounterSet() {
		LocusCategory[] locusCategories = LocusCategory.values();
		categoryCount = locusCategories.length;
		categoriesByValue = new HashMap<LocusCategory,IntegerCounter>(categoryCount);
		categories = new IntegerCounter[locusCategories.length];
		int nextIndex = 0;
		for (LocusCategory lc : locusCategories)
			categoriesByValue.put(lc, categories[nextIndex++] = new IntegerCounter());
		all = new IntegerCounter();
	}
	
	IntegerCounter getCategoryCounter(LocusCategory c) {
		return categoriesByValue.get(c);
	}
	
	void addValue(int value, byte categories) {
		all.addValue(value);
		byte mask = 1;
		for (int i = 0; i < categoryCount; i++) {
			if ((categories & mask) != 0) 
				this.categories[i].addValue(value);
			mask <<= 1;
		}		
	}
	
	@Deprecated
	void applyIncrement(int depth, byte categories) {
		addValue(depth,categories);
	}
	
	public void applyCounterSet(IntegerSequenceCounterSet rhs) {
		all.applyCounter(rhs.all);
		for (int i = 0; i < categoryCount; i++) 
			categories[i].applyCounter(rhs.categories[i]);
	}

	public void printReport(String indent, PrintWriter pw) {
		all.printReport(indent + "  ",pw);
		for (Map.Entry<LocusCategory,IntegerCounter> e : categoriesByValue.entrySet()) {
			if (e.getValue().toDistribution().count() == 0) continue;
			pw.println(indent + "Locus Category " + e.getKey() + ":");
			e.getValue().printReport(indent + "  ",pw);
		}
	}
	
	

}
