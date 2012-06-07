package net.malariagen.gatk.genotyper;

public enum GenotypeVariantFilterEmitMode {
	/** 
	 * Do not filter variants, show them all. 
	 */
	EMIT_ALL,
	/**
	 * Filter does variant where there is no genotype different than than the reference.
	 */
	VARIANT, 
	/**
	 * Filter does variants where there are no different genotypes between samples (even if all are non-reference).
	 */
	POLYMORPHIC, 
	/**
	 * Emit all variants but indicate lack of variant genotypes using a filter. 
	 */
	NO_VARIANT_FILTER,
	/**
	 * Emit all variants but indicate lack of polymorphism with a filter. 
	 */
	NO_POLYMORPHIC_FILTER;
}
