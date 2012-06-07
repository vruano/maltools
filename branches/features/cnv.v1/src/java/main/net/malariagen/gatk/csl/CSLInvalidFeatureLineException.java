package net.malariagen.gatk.csl;

public class CSLInvalidFeatureLineException extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public CSLInvalidFeatureLineException(String line, String message) {
		super((message == null ? "" : message)  + ": " + line);
	}

}
