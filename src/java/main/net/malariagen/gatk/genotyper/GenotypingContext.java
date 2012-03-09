package net.malariagen.gatk.genotyper;

import java.util.List;

import org.broadinstitute.sting.utils.variantcontext.Allele;

/**
 * Collects general information about the context in where genotyping is going to be performed.
 */
public interface GenotypingContext {

	/**
	 * Returns the reference allele in this context.
	 * 
	 * @return never {@code null} but perhaps a {@link Allele#NO_CALL} allele reference to indicate the absence of a reference allele.
	 */
	public Allele getReferenceAllele();
	
	/**
	 * Number of possible alleles in this context. At any rate the reference allele (perhaps a {@link Allele#NO_CALL}) is always counted.
	 * @return one or greater.
	 */
	public int getAlleleCount();
	
	/**
	 * Returns the i-th plausible allele at this context. By convention the reference is has index 0.  
	 * @throws IllegalArgumentException if {@code index} is less than 0 or greater than {@code {@link #getAlleleCount()} - 1}.
	 */
	public Allele getAllele(int index);
	
	/**
	 * Returns the plausible allele based on the base byte.
	 * 
	 * @return {@code null} if there is no allele with the base indicate with {@code b}. 
	 */
	public Allele getAllele(byte b);

	/**
	 * Return the index of the allele give its byte.
	 * @param b
	 * @return -1 if there is not such allele.
	 */
	public int getAlleleIndex(byte b);
	
	/**
	 * List of all possible alleles sorted by their index. The first allele is always the reference. 
	 * 
	 * @throws never {@code null} nor empty (the reference allele is always present).
	 */
	public List<Allele> getAlleleList();
		
	public boolean hasAllele(Allele a);
	
	public boolean hasAllele(byte b);
}
