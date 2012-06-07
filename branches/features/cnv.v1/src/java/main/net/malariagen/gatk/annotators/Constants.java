package net.malariagen.gatk.annotators;

import java.text.NumberFormat;

public class Constants {
	
	public static final String MINOR_READ_ALLELE_FREQUENCY_KEY = "MrAF";

	public static final String MINOR_READ_ALLELE_MAXIMUM_SAMPLE_DEPTH_KEY = "MrAmaxSD";

	public static final String MINOR_READ_ALLELE_TOTAL_DEPTH_KEY = "MrAD";

	public static final String MINOR_ALLELE_KEY = "MA";
	
	public static final String MINOR_ALLELE_FREQUENCY_KEY = "MAF";	
	
	public static final String MINOR_READ_ALLELE_KEY = "MrA";

	public static final String ABSOLUTE_ALLELE_MAX_SAMPLE_DEPTH_KEY = "AmaxSD";

	public static final String ABSOLUTE_TOTAL_DEPTH_KEY = "AT";

	public static final String ABSOLUTE_ALLELE_DEPTH_KEY = "AD";

	public static final String UNIQUENESS_KEY = "UQ";

	public static final String UNIQUENESS_ROD_NAME = "uniqueness";

	public static final String FEATURES_ROD_NAME = "features";

	public static final String CODING_KEY = "CODING";
	
	public static final String LOWEST_MAXIMUM_READ_SAMPLE_DEPTH_KEY = "LMSrSD";

	public static final String COVERAGE_MEDIAN_FOLD_KEY = "CMF";
	
	public static final String COVERAGE_CUMULATIVE_PROBABILITY_KEY = "CCP";
	
	public static final int FRACTION_ANNOTATION_PRECISION = 4;
	
	private static final NumberFormat FRACTION_ANNOTATION_FORMATTER;
	
	public static final String READ_COUNTS_KEY = "RC";
	
	
	static {
		NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(FRACTION_ANNOTATION_PRECISION);
		FRACTION_ANNOTATION_FORMATTER = formatter;
	}
	
	public static String formatFractionAnnotation(double f) {
		return FRACTION_ANNOTATION_FORMATTER.format(f);
	}

	
}
