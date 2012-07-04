package net.malariagen.gatk.coverage;

import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

import net.malariagen.gatk.coverage.SequenceComplexity.LocusComplexity;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;

public class MultiWindowSequenceComplexity {

	protected SequenceComplexity[] complexities;
	protected Map<Integer, SequenceComplexity> byWs;
	protected int maxWindow = 0;

	protected PriorityQueue<WindowSet> windows = new PriorityQueue<WindowSet>();
	protected Map<GenomeLoc, WindowSet> windowByLoc = new HashMap<GenomeLoc, WindowSet>(
			100);

	public List<Map<Integer, LocusComplexity>> count(ReferenceContext ref) {

		GenomeLoc loc = ref.getLocus();
		WindowSet ws = windowByLoc.get(loc);
		if (ws == null) {
			if (!windows.isEmpty() && !windows.element().start.getContig().equals(loc.getContig())) {
				windowByLoc.clear();
				windows.clear();
			}
			windowByLoc.put(loc, ws = new WindowSet(loc));
			windows.add(ws);
		}

		for (Integer i : byWs.keySet()) {
			LocusComplexity lc = byWs.get(i).count(ref);
			if (lc == null)
				continue;
			ws = windowByLoc.get(lc.getLocus());
			if (ws == null)
				throw new RuntimeException("complexity at locus " + loc
						+ " seen before locus being visited!!!");
			ws.bySize.put(i, lc);
		}
		if (windows.element().bySize.size() == byWs.size()) {
			List<Map<Integer, LocusComplexity>> result = new LinkedList<Map<Integer, LocusComplexity>>();
			while (!windows.isEmpty()) {
				WindowSet ws2 = windows.element();
				if (ws2.bySize.size() < byWs.size())
					break;
				result.add(ws2.bySize);
				windowByLoc.remove(ws2.start);
				windows.remove();
			}
			return result;
		}
		else {
		  return Collections.emptyList();
		}
	}

	MultiWindowSequenceComplexity(int... ws) {
		for (int i : ws)
			if (i <= 0)
				throw new IllegalArgumentException(
						"the requiested windows size cannot be less or equal to 0");
		byWs = new HashMap<Integer, SequenceComplexity>(ws.length);
		for (int i : ws) {
			if (i > maxWindow)
				maxWindow = i;
			if (byWs.get(i) == null)
				byWs.put(i, SequenceComplexity.create(i));
		}

		complexities = byWs.values().toArray(
				new SequenceComplexity[byWs.size()]);
	}

	protected class WindowSet implements Comparable<WindowSet> {
		GenomeLoc start;
		Map<Integer, LocusComplexity> bySize;

		public WindowSet(GenomeLoc loc) {
			start = loc;
			bySize = new HashMap<Integer, LocusComplexity>();
		}

		public boolean equal(Object o) {
			if (o == null)
				return false;
			else if (o instanceof WindowSet)
				return equal((WindowSet)o);
			else 
				return false;
		}
		
		public int hashCode() {
			return start.hashCode();
		}
		
		public boolean equal(WindowSet o) {
			return compareTo(o) == 0;
		}

		@Override
		public int compareTo(WindowSet o) {
			return start.compareTo(o.start);
		}

	}

}
