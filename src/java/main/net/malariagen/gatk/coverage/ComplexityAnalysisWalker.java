package net.malariagen.gatk.coverage;

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.coverage.SequenceComplexity.LocusComplexity;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

public class ComplexityAnalysisWalker extends LocusWalker<ReferenceContext,SequenceComplexity> {

	@Argument(shortName="W", fullName="windowSize", required=true, doc = "window-size to get stats on")
	private int windowSize = 0;
	
	@Output(shortName = "o", fullName = "output", doc = "File to which variants should be written", required = true)
	protected VCFWriter writer = null;

	@Override
	public ReferenceContext map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		return ref;
	}
	
	@Override
	public void initialize() {
		super.initialize();
		if (writer == null)
			throw new IllegalStateException("the write is yet null");
		if (windowSize <= 0)
			throw new IllegalArgumentException("a strictly positive window size must be indicated, the one provided was " + windowSize + ")");
		
		Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		headerLines.add(new VCFInfoHeaderLine("GCBias", 1, VCFHeaderLineType.Float, "GC bias expressed in a percentage from 0 to 100"));
		headerLines.add(new VCFInfoHeaderLine("NucEnt",1, VCFHeaderLineType.Float,"Nucleotide entropy using Euler's e as base"));
		headerLines.add(new VCFInfoHeaderLine("GCHet",1,VCFHeaderLineType.Float,"GC bias heterogeneity"));
		headerLines.add(new VCFInfoHeaderLine("END",1,VCFHeaderLineType.Float,"Stop positio of the interval considerde in complex INFO fields"));
		VCFHeader header = new VCFHeader(headerLines);
		writer.writeHeader(header);
	}
	
	@Override
	public SequenceComplexity reduceInit() {
		return SequenceComplexity.create(windowSize);
	}

	@Override
	public SequenceComplexity reduce(ReferenceContext value, SequenceComplexity sum) {
		SequenceComplexity.LocusComplexity lc = sum.count(value);
		if (lc != null)
			emit(lc);
		return sum;
	}

	private void emit(LocusComplexity lc) {
		VariantContextBuilder vcb = new VariantContextBuilder();
		Collection<Allele> alleles = Collections.singleton(Allele.create(lc.getRefNuc().byteValue(),true));
		vcb.alleles(alleles);
		Map<String,Object> attributes = new LinkedHashMap<String,Object>(4);
		attributes.put("GCBias",lc.getGcBias());
		attributes.put("NucEnt",lc.getNucEnt());
		attributes.put("GCHet",lc.getGcHet());
		attributes.put("END",lc.getLocus().getStart() + windowSize - 1);
		vcb.loc(lc.getLocus());
		vcb.attributes(attributes);
		VariantContext vc = vcb.make();

		writer.add(vc);
	}

}
