package net.malariagen.gatk.gff;

public class GFFInvalidFeatureLineException extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 5000119694654711590L;

	public GFFInvalidFeatureLineException(String line, String message) {
		super((message == null ? "" : message)  + ": " + line);
	}

}
