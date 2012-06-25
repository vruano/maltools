package net.malariagen.gatk.genotyper;

import java.util.Collections;
import java.util.List;

import net.sf.picard.cmdline.CommandLineParser;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;

public class MetaArgumentCollection extends UnifiedArgumentCollection {

    // control the various models to be used
    @Argument(fullName = "speciment_model", shortName = "smodel", doc = "Per specifiment probabilistic model", required = false)
    public String smodel = "Bernoulli(0.5)";
    
    @Argument(fullName = "baseq_dist_output", shortName = "baseq_do", doc = "Base quality distribution output file", required = false)
    public String baseqDistOut = null;

    @Argument(fullName = "min_genotype_quality", shortName = "mgq", doc = "Minimum genotype quality to make a call", required = false)
	public double MIN_GENOTYPE_QUALITY = - 0.01;
    
    @Argument(fullName = "gt_emit_mode", shortName = "gem", doc = "Wether to emit variant based on their genotype calls", required = false)
    public GenotypeVariantFilterEmitMode gtVarFilterEmitMode = GenotypeVariantFilterEmitMode.EMIT_ALL;
    
    @Argument(fullName = "min_genotype_confidence", shortName = "mgc", doc = "Minimum genotype confidence to make a call", required = false)
    public double MIN_GENOTYPE_CONFIDENCE = -100;
    
    @Argument(fullName = "parent", shortName = "P", doc = "Parent sample name", required = false)
    public List<String> parents = Collections.emptyList();
   
    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum mapping  quality required for a read to be considered during genotyping", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 0;
    
    @Argument(fullName = "rodBind", shortName = "B", doc = "Additional ROD bindings", required = false)
    public List<RodBinding<Feature>> rodBinds;
    
}
