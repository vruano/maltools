package net.malariagen.gatk.walker;

public enum FragmentFilter {

	NONE, GENOMIC_MASK, INVERSION, TRANSLOCATION, FIRST_FACE_AWAY, LAST_FACE_AWAY, BOTH_FACE_AWAY, UNMAPPED, BOTH_UNMAPPED, FIRST_UNMAPPED, LAST_UNMAPPED,
	BOTH_LOW_MQ, FIRST_LOW_MQ, LAST_LOW_MQ, LOW_MQ, INSERTION, DELETION, PROPER_PAIR_FLAG, STRUCTURAL_VARIATION;
	
	
	private static int[] IMPLIES_MASK;
	
	static {
		
		IMPLIES_MASK = new int[values().length];

		for (int i = 0; i < values().length; i++) 
			IMPLIES_MASK[i] = 1 << i;
		
		IMPLIES_MASK[BOTH_FACE_AWAY.ordinal()] |=  impliesMask(BOTH_FACE_AWAY,INVERSION,FIRST_FACE_AWAY,LAST_FACE_AWAY,STRUCTURAL_VARIATION);
		IMPLIES_MASK[FIRST_FACE_AWAY.ordinal()] |= impliesMask(FIRST_FACE_AWAY,INVERSION,STRUCTURAL_VARIATION);
		IMPLIES_MASK[LAST_FACE_AWAY.ordinal()] |= impliesMask(LAST_FACE_AWAY,INVERSION,STRUCTURAL_VARIATION);

		IMPLIES_MASK[BOTH_LOW_MQ.ordinal()] |=  impliesMask(BOTH_LOW_MQ,LOW_MQ,FIRST_LOW_MQ,LAST_LOW_MQ);
		IMPLIES_MASK[FIRST_LOW_MQ.ordinal()] |= impliesMask(FIRST_LOW_MQ,LOW_MQ);
		IMPLIES_MASK[LAST_LOW_MQ.ordinal()] |= impliesMask(LAST_LOW_MQ,LOW_MQ);

		IMPLIES_MASK[BOTH_UNMAPPED.ordinal()] |=  impliesMask(BOTH_UNMAPPED,UNMAPPED,FIRST_UNMAPPED,LAST_UNMAPPED) | IMPLIES_MASK[BOTH_LOW_MQ.ordinal()];
		IMPLIES_MASK[FIRST_UNMAPPED.ordinal()] |= impliesMask(FIRST_UNMAPPED,UNMAPPED) | IMPLIES_MASK[FIRST_LOW_MQ.ordinal()];
		IMPLIES_MASK[LAST_UNMAPPED.ordinal()] |= impliesMask(LAST_UNMAPPED,UNMAPPED) | IMPLIES_MASK[LAST_LOW_MQ.ordinal()];
	}
	
	private static int impliesMask(FragmentFilter ... filters) {
		int result = 0;
		for (FragmentFilter filter : filters)
			result |= 1 << filter.ordinal();
		return result;
	}
	
	public boolean implies(FragmentFilter other) {
		if (other == this)
			return true;
		if (other == null)
			throw new IllegalArgumentException("other cannot be null");
		int impliesMask = IMPLIES_MASK[this.ordinal()];
		
		return (impliesMask & (1 << other.ordinal())) != 0;
	}
	
}
