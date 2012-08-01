package net.malariagen.gatk.utils;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.broadinstitute.sting.commandline.Argument;

import net.malariagen.gatk.walker.FragmentLengthSummary;

public class WindowSizeDB {
	
	private Integer defaultWindowSize;
	private Map<String, Integer> perSampleOrGroup;
	
	@SuppressWarnings("unchecked")
	WindowSizeDB(Integer defaultWindowSize, Map<String,Integer> perSampleOrGroup) {
		this.defaultWindowSize = defaultWindowSize;
		this.perSampleOrGroup = perSampleOrGroup == null ? Collections.EMPTY_MAP : perSampleOrGroup;
	}
	
	public int getWindowSize(String readGroupOrSample) {
		if (readGroupOrSample == null)
			return getDefaultWindowSize();
		else if (!perSampleOrGroup.containsKey(readGroupOrSample)) 
			return getDefaultWindowSize();
		Integer result = perSampleOrGroup.get(readGroupOrSample);
		return result == null ? getDefaultWindowSize() : result;
	}

	public int getDefaultWindowSize() {
		return defaultWindowSize == null ? -1 : defaultWindowSize;
	}
	
	public static WindowSizeDB create(Integer defaultWindowSize, FragmentLengthSummary fls) {
		Map<String,Integer> perSampleOrGroup = new HashMap<String,Integer>(fls.getReadGroups().size() + fls.getSamples().size());
		for (String rg : fls.getReadGroups()) 
			perSampleOrGroup.put(rg,(int) Math.round(fls.getSampleFragmentLengths(rg).median()));
		
		for (String sm : fls.getSamples()) 
			perSampleOrGroup.put(sm,(int) Math.round(fls.getSampleFragmentLengths(sm).median()));
		return new WindowSizeDB(defaultWindowSize,perSampleOrGroup);
	}
	
}
