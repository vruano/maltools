package net.malariagen.gatk.cnv;

import org.broadinstitute.sting.utils.GenomeLoc;

public interface CopyNumberStatistics {

	public GenomeLoc getLocation();

	/**
	 * Returns the read-starts count at #offset. 
	 * 
	 * @param offset
	 * @return
	 */
	public int getReadStarts(int offset); 
	
	/**
	 * @returns reference to the array containing the read start count of the
	 *          region. Please do not change its content.
	 */
	public int[] getReadStarts();

}
