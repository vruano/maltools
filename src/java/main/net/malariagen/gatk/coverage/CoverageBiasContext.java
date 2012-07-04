package net.malariagen.gatk.coverage;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

class CoverageBiasContext extends VariantContext {

    private ReferenceContext refContext;	
	
	protected CoverageBiasContext(VariantContext other, ReferenceContext ref) {
		super(other);
		refContext = ref;
	}
	
	ReferenceContext getReferenceContext() {
		return refContext;
	}

}
