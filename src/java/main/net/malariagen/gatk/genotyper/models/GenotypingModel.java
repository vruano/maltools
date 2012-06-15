package net.malariagen.gatk.genotyper.models;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import net.malariagen.gatk.genotyper.GenotypingContext;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.MutableGenotype;

public interface GenotypingModel {

	/**
	 * Return the genotying model name.
	 * @return never {@code null}
	 */
	public String getModelName();
	
	/**
	 * Indicates how many possible genotypes are there for a speciment or sample.
	 * @return zero or greater however it should typically have at least two possible genotypes for the model to be informative.
	 */
	public int getGenotypeCount();
	
	/**
	 * Returns the set of prior for each possible genotype.
	 */
	public GenotypePriors getGenotypePriors();
	
	public void setGenotypingContext(GenotypingContext gc);
	
	/**
	 * Returns the genotyping context
	 * @return never {@code null}
	 */
	public GenotypingContext getGenotypingContext();

	/**
	 * Returns the list of alleles for the i-th genotype.
	 * @param 
	 * @param genotype must be a valid genotype index for this model, thus between {@code 0} and {@code {@link #getGenotypeCount()}-1}
	 * @return never {@code null} but possibly a zero-length list meaning the trivial no-call genotype.
	 */
	public List<Allele> getGenotypeAlleles(int genotype);

	public Map<String, MutableGenotype> callGenotypes(Map<String, AlignmentContext> sc);
	
	/**
	 * Returns the variant quality given the stratified alignment content, genotypes and posteriors.
	 * @param ac
	 * @return
	 */
	public double variantQuality(Map<String, AlignmentContext> ac, Map<String, Genotype> gt);

	/**
	 * Returns the genotype index from the allele list provided. 
	 * 
	 * @param alleles
	 * @return -1 if this list is not possible for this model.
	 */
	public int getGenotypeIndex(List<Allele> alleles);


	public Collection<? extends VCFHeaderLine> getHeaderLines();
	
}
