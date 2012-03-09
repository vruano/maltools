package net.malariagen.gatk.genotyper;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;

public class MetaArgumentCollection extends UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "speciment_model", shortName = "smodel", doc = "Per specifiment probabilistic model", required = false)
    public String smodel = "Bernoulli(0.5)";
    
    @Argument(fullName = "baseq_dist_output", shortName = "baseq_do", doc = "Base quality distribution output file", required = false)
    public String baseqDistOut = null;

    @Argument(fullName = "min_genotype_quality", shortName = "mgq", doc = "Minimum genotype quality to make a call", required = false)
	public double MIN_GENOTYPE_QUALITY = 0.01;
    
    @Argument(fullName = "gt_emit_mode", shortName = "gem", doc = "Wether to emit variant based on their genotype calls", required = false)
    public GenotypeVariantFilterEmitMode gtVarFilterEmitMode = GenotypeVariantFilterEmitMode.EMIT_ALL;
    
    @Argument(fullName = "min_genotype_confidence", shortName = "mgc", doc = "Minimum genotype confidence to make a call", required = false)
    public double MIN_GENOTYPE_CONFIDENCE = -100;
    
}
