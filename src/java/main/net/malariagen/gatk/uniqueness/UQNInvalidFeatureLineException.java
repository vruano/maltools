package net.malariagen.gatk.uniqueness;

public class UQNInvalidFeatureLineException extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public UQNInvalidFeatureLineException(String line, String message) {
		super((message == null ? "" : message)  + ": " + line);
	}

}
