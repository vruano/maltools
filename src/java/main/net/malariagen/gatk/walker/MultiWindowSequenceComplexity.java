package net.malariagen.gatk.walker;

import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;


import net.malariagen.gatk.walker.ReferenceComplexityWalker.Whence;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;

public class MultiWindowSequenceComplexity {

	protected SequenceComplexity[] complexities;
	protected Map<Integer, SequenceComplexity> byWs;
	protected int maxWindow = 0;
	protected Whence whence;

	protected PriorityQueue<WindowSet> windows = new PriorityQueue<WindowSet>();
	protected Map<GenomeLoc, WindowSet> windowByLoc = new HashMap<GenomeLoc, WindowSet>(
			100);
	private GenomeLoc lastLocus;
        private GenomeLoc lastEmited;
	
	public List<Map<Integer, SequenceComplexity.LocusComplexity>> flush() {
		for (Integer i : byWs.keySet()) 
			for (SequenceComplexity.LocusComplexity slc : byWs.get(i).flush()) {
				GenomeLoc loc = slc.getLocus();
                                if (windowByLoc.get(loc) == null)
                                   System.err.println(" " + lastLocus + " " + lastEmited + " " + loc);
				windowByLoc.get(loc).bySize.put(i,slc);
                        }
		List<Map<Integer, SequenceComplexity.LocusComplexity>> result = new LinkedList<Map<Integer, SequenceComplexity.LocusComplexity>>();
		while (!windows.isEmpty()) {
			WindowSet ws2 = windows.remove();
			if (ws2.bySize.size() < byWs.size())
				throw new IllegalStateException("some window size did not flush to the end");
			result.add(ws2.bySize);
			windowByLoc.remove(ws2.start);
		}
		return result;
		
	}
	
	public List<Map<Integer, SequenceComplexity.LocusComplexity>> count(ReferenceContext ref) {

		GenomeLoc loc = ref.getLocus();
		lastLocus = loc;
		WindowSet ws = windowByLoc.get(loc);
		if (ws == null) {
			if (!windows.isEmpty() && !windows.element().start.getContig().equals(loc.getContig())) {
		        throw new IllegalArgumentException("what tha");
//				windows.clear();
			}
			windowByLoc.put(loc, ws = new WindowSet(loc));
			windows.add(ws);
		}
		
		for (Integer i : byWs.keySet()) {
			SequenceComplexity.LocusComplexity lc = byWs.get(i).count(ref);
			if (lc == null)
				continue;
			ws = windowByLoc.get(lc.getLocus());
			if (ws == null)
				throw new RuntimeException("complexity at locus " + loc
						+ " seen before locus being visited!!!");
			ws.bySize.put(i, lc);
		}
		if (windows.element().bySize.size() == byWs.size()) {
			List<Map<Integer, SequenceComplexity.LocusComplexity>> result = new LinkedList<Map<Integer, SequenceComplexity.LocusComplexity>>();
			while (!windows.isEmpty()) {
				WindowSet ws2 = windows.element();
				if (ws2.bySize.size() < byWs.size())
					break;
				result.add(ws2.bySize);
				windowByLoc.remove(ws2.start);
				windows.remove();
                                lastEmited = ws2.start;
			}
			return result;
		}
		else {
		  return Collections.emptyList();
		}
	}

	MultiWindowSequenceComplexity(int[] ws, Whence whence) {
		this.whence = whence;
		for (int i : ws)
			if (i <= 0)
				throw new IllegalArgumentException(
						"the requiested windows size cannot be less or equal to 0");
		byWs = new HashMap<Integer, SequenceComplexity>(ws.length);
		for (int i : ws) {
			if (i > maxWindow)
				maxWindow = i;
			int padding;
			switch (whence) {
			case START: padding = 0; break;
			case END: padding = i - 1; break;
			default: // CENTER
				padding = i >> 1;
			}
			if (byWs.get(i) == null) 
				byWs.put(i, SequenceComplexity.create(i,padding));
		}

		complexities = byWs.values().toArray(
				new SequenceComplexity[byWs.size()]);
	}

	protected class WindowSet implements Comparable<WindowSet> {
		GenomeLoc start;
		Map<Integer, SequenceComplexity.LocusComplexity> bySize;

		public WindowSet(GenomeLoc loc) {
			start = loc;
			bySize = new HashMap<Integer, SequenceComplexity.LocusComplexity>();
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

	public GenomeLoc lastLocus() {
		return lastLocus;
	}

}
