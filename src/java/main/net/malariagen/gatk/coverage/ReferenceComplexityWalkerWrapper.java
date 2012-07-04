package net.malariagen.gatk.coverage;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

public class ReferenceComplexityWalkerWrapper {

	private GenomeAnalysisEngine toolkit;
	
	public ReferenceComplexityWalkerWrapper(GenomeAnalysisEngine gae) {
		this.toolkit = gae;
	}
	
	
	
	public void initialize() {
		
	}
}
