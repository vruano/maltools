package net.malariagen.gatk.coverage;

class AccessibilityStatistics {

	public AccesibilityPerLocusCategory CODING = new AccesibilityPerLocusCategory();
	public AccesibilityPerLocusCategory NON_CODING = new AccesibilityPerLocusCategory();
	
	
	class AccesibilityPerLocusCategory {
		public long total;
	
		public long accessibleByUniqueness;
	
		public long accessibleByCoverage;
	
		public long accessibleByAny;
	
		public long nonAccessibleByUniqueness;
	
		public long nonAccessibleDueToLowCoverage;
	
		public long nonAccessibleDueToHighCoverage;
	}
}
