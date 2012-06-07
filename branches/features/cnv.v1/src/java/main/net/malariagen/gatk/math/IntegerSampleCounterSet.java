package net.malariagen.gatk.math;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import net.malariagen.gatk.coverage.LocusCategory;


class IntegerSampleCounterSet {

	final Map<String,IntegerSequenceCounterSet> sequencesByName;
	final IntegerSequenceCounterSet[] sequences;
	final IntegerSequenceCounterSet all;
	private final int sequenceCount;
	
	IntegerSampleCounterSet(String ... sequenceNames) {
		sequenceCount = sequenceNames.length;
		sequencesByName = new HashMap<String,IntegerSequenceCounterSet>(sequenceCount);
		sequences = new IntegerSequenceCounterSet[sequenceCount];
		int nextIndex = 0;
		for (String seq : sequenceNames)
			sequencesByName.put(seq, sequences[nextIndex++] = new IntegerSequenceCounterSet());
		all = new IntegerSequenceCounterSet();
	}
	
	IntegerSequenceCounterSet getCategoryCounter(LocusCategory c) {
		return sequencesByName.get(c);
	}
	
	void addValue(int value, byte categories, int sequence) {
		all.addValue(value,categories);
		sequences[sequence].addValue(value,categories);
	}
	
	@Deprecated
	void applyIncrement(int depth, byte categories, String sequence) {
		all.addValue(depth,categories);
		sequencesByName.get(sequence).addValue(depth,categories);
	}

	public void applyCounterSet(IntegerSampleCounterSet rhs) {
		all.applyCounterSet(rhs.all);
		for (int i = 0; i < sequenceCount; i++) 
			sequences[i].applyCounterSet(rhs.sequences[i]);
	}

	public void printReport(String indent, PrintWriter pw) {
		pw.print(indent + "Coverage across all sequences:");
		all.printReport(indent + "  ",pw);
		for (Map.Entry<String,IntegerSequenceCounterSet> e : sequencesByName.entrySet() ) {
			if (e.getValue().all.toDistribution().count() == 0) continue;
			pw.println(indent + "Sequence " + e.getKey() + ":");
			e.getValue().printReport(indent + "  ",pw);
		}
	}

	

}
