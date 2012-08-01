package net.malariagen.gatk.walker;

public enum FragmentOrder {

	/**
	 * Fragment first read location in the genome. 
	 */
	FRAGMENT_START, 
	
	/**
	 * First fragment position in the reference genome.
	 */
	REFERENCE_START, 
	
	/**
	 * No sorting required.
	 */
	NONE;
	
	
}
