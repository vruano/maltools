package net.malariagen.gatk.walker;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.coverage.CoverageBiasWalker.GroupBy;
import net.malariagen.gatk.walker.SequenceComplexity.LocusComplexity;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

@By(DataSource.REFERENCE)
@PartitionBy(PartitionType.CONTIG)
@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
public class ReferenceComplexityWalker
		extends
		LocusWalker<net.malariagen.gatk.walker.ReferenceComplexityWalker.Locus, MultiWindowSequenceComplexity> {

	public enum Whence {
		START,CENTER,END;
	}

	static class Locus {
		ReferenceContext ref;
	}

	@Argument(shortName = "W", fullName = "windowSize", required = false, doc = "window-size to get stats on")
	protected List<Integer> windowSizeList = Collections.emptyList();

	@Argument(shortName = "useInfoFields", fullName = "useInfoFields", required = false, doc = "whether complexity for a fixed unique window size shoul be outputted in the info fields")
	protected boolean useInfoFields = false;

	protected Integer windowSize;

	@Argument(shortName = "rounding", fullName = "windowSizeRounding", doc = "window size rounding, 1 means no rounding", required = false)
	protected int rounding = 1;

	@Output(shortName = "o", fullName = "output", doc = "File to which variants should be written", required = true)
	protected VCFWriter writer = null;
	
	@Argument(shortName = "whence", fullName="whence", doc ="Where the targeted position in the reference is located within the window (START,CENTER or END) ", required = false)
	protected Whence whence = Whence.CENTER;

	@Override
	public Locus map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		Locus result = new Locus();
		result.ref = ref;
		return result;
	}

	@Output(shortName = "groupBy", required = false)
	protected GroupBy groupBy = null;

	Map<String, Integer> groupWindowSize;
	Map<String, Integer> realGroupWindowSize;


	@Override
	public void initialize() {
		super.initialize();
		if (writer == null)
			throw new IllegalStateException("the write is yet null");
		if (!windowSizeList.isEmpty())
			windowSize = windowSizeList.get(0);
		else
			throw new IllegalArgumentException("you must indicate at least one window size (-W)");
		if (rounding <= 0)
			throw new UserException("the window size rounding must be positive");

		if (useInfoFields && windowSizeList.size() != 1)
			throw new UserException(
					"you cannot request complexity to be output in info fields unless you indicate one and only one window size with -W");

		Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		VCFHeader header = null;

		groupWindowSize = new HashMap<String, Integer>(100);
		for (Integer ws : this.windowSizeList) {
			if (ws == null || ws <= 0)
				throw new UserException("no requested window size can be null, zero or negative (-W " + ws + ")");
			int actualWs = (int) Math.round(ws/rounding) * rounding;
			groupWindowSize.put("" + ws, actualWs);
		}
		Set<String> groupNames;
		if (useInfoFields) {
			groupNames = Collections.emptySet();
			headerLines.add(new VCFInfoHeaderLine("GCBias", 1,
					VCFHeaderLineType.Float,
					"GC bias expressed in a percentage from 0 to 100 (only if NucCnt > 0)"));
			headerLines.add(new VCFInfoHeaderLine("NucEnt", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy in nats (only if NucCnt > 0)"));
			headerLines.add(new VCFInfoHeaderLine("TriEnt", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy in nats (only if TriCnt > 0)"));
			headerLines.add(new VCFInfoHeaderLine("NucCnt", 1,
					VCFHeaderLineType.Integer,
					"Number of sense nucleotides in the window (A, C, G or T)"));
			headerLines.add(new VCFInfoHeaderLine("TriCnt", 1,
					VCFHeaderLineType.Integer,
					"Number of sense trinucleotides in the window (composed only by A, C, G or T)"));
			headerLines.add(new VCFInfoHeaderLine("DUST", 1,
					VCFHeaderLineType.Float,
					"DUST complexity metric for the given window (only if TriCnt > 1)"));
			switch (whence) {
			case CENTER:
				headerLines.add(new VCFInfoHeaderLine("START", 1,
						VCFHeaderLineType.Float,
						"Stat position of the window interval"));
			case START:
				headerLines.add(new VCFInfoHeaderLine("END", 1,
					VCFHeaderLineType.Float,
					"Stop position of the window interval"));
				break;
			case END:
				headerLines.add(new VCFInfoHeaderLine("START", 1,
						VCFHeaderLineType.Float,
						"Start position of the window interval"));
			}
		
		}
		else {
			groupNames = groupWindowSize.keySet();
			headerLines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String,
					"Genotype Call, never used but left in to be compliant with all standards"));
			headerLines.add(new VCFFormatHeaderLine("GC", 1,
					VCFHeaderLineType.Float,
					"GC bias expressed in a percentage from 0 to 100 (only if NC > 0)"));
			headerLines.add(new VCFFormatHeaderLine("NE", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy in nats (only if NC > 0)"));
			headerLines.add(new VCFFormatHeaderLine("TE", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy in nats (only if TC > 0)"));
			headerLines.add(new VCFFormatHeaderLine("NC", 1,
					VCFHeaderLineType.Integer,
					"Number of sense nucleotides in the window (A, C, G or T)"));
			headerLines.add(new VCFFormatHeaderLine("TC", 1,
					VCFHeaderLineType.Integer,
					"Number of sense trinucleotides in the window (composed only by A, C, G or T)"));
			headerLines.add(new VCFFormatHeaderLine("DU", 1,
					VCFHeaderLineType.Float,
					"DUST complexity metric for the given window (only if TC > 1)"));
			headerLines.add(new VCFFormatHeaderLine("ED", 1,
					VCFHeaderLineType.Float,
					"Stop position of the window interval"));
			switch (whence) {
			case CENTER:
				headerLines.add(new VCFInfoHeaderLine("ST", 1,
						VCFHeaderLineType.Float,
						"Start position of the window interval"));
			case START:
				headerLines.add(new VCFInfoHeaderLine("ED", 1,
					VCFHeaderLineType.Float,
					"Stop position of the window interval"));
				break;
			case END:
				headerLines.add(new VCFInfoHeaderLine("ST", 1,
						VCFHeaderLineType.Float,
						"Start position of the window interval"));
			}
			
		}
		header = new VCFHeader(headerLines,groupNames);
		writer.writeHeader(header);
	}


	@Override
	public MultiWindowSequenceComplexity reduceInit() {
		Set<Integer> windowSizes = new HashSet<Integer>();
		if (groupBy != GroupBy.NONE)
			windowSizes.addAll(groupWindowSize.values());
		if (windowSize != null)
			windowSizes.add(windowSize);
		int[] windowSizeInts = new int[windowSizes.size()];
		int nextIdx = 0;
		for (Integer i : windowSizes) {
			windowSizeInts[nextIdx++] = i;
		}
		return new MultiWindowSequenceComplexity(windowSizeInts,whence);
	}

	@Override
	public MultiWindowSequenceComplexity reduce(Locus value,
			MultiWindowSequenceComplexity sum) {
		if (sum.lastLocus() != null && value.ref.getLocus().compareContigs(sum.lastLocus()) != 0) 
			for (Map<Integer, SequenceComplexity.LocusComplexity> l : sum.flush())
				if (l.size() != 0) emit(l);
		sum = this.reduceInit();
		List<Map<Integer, SequenceComplexity.LocusComplexity>> lcm = sum.count(
				value.ref);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : lcm)
			if (l.size() != 0)
				emit(l);
		return sum;
	}
	
	@Override
	public void onTraversalDone(MultiWindowSequenceComplexity sum) {
		super.onTraversalDone(sum);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : sum.flush())
			if (l.size() != 0) emit(l);
	}

	private void emit(Map<Integer, LocusComplexity> lcm) {
		VariantContextBuilder vcb = new VariantContextBuilder();
		LocusComplexity example = lcm.values().iterator().next();
		Collection<Allele> alleles = Collections.singleton(Allele.create(
				example.getRefNuc().byteValue(), true));
		vcb.alleles(alleles);
		GenomeLoc loc = example.getLocus();
		Map<String, Object> attributes = new LinkedHashMap<String, Object>(4);
		if (useInfoFields && lcm.containsKey(windowSize)) {
			LocusComplexity lc = lcm.get(windowSize);
			if (lc.getNucCount() > 0) {
				attributes.put("NucEnt", lc.getNucEnt());
				attributes.put("GCBias", lc.getGcBias());
			}
			if (lc.getTrinucCount() > 0) {
				attributes.put("TriEnt", lc.getTriEnt());
				if (lc.getTrinucCount() > 1) attributes.put("DUST", lc.getDUST());
			}
			attributes.put("NucCnt", lc.getNucCount());
			attributes.put("TriCnt", lc.getTrinucCount());
			switch (whence) {
			case CENTER:
				attributes.put("START",lc.getLocus().getStart());
			case START:  attributes.put("END", lc.getLocus().getStart() + lc.size() - 1); break;
			case END:
				attributes.put("START",lc.getLocus().getStart());
			}
			if (lc.getLocus().compareTo(loc) < 0)
				loc = lc.getLocus();
		}
		vcb.attributes(attributes);
		List<Allele> noCall = Collections.singletonList(Allele.NO_CALL);
		Set<String> noFilters = Collections.emptySet();
		GenotypesContext gc = GenotypesContext.create();
		for (Map.Entry<String, Integer> gws : groupWindowSize.entrySet()) {
			LocusComplexity lc = lcm.get(gws.getValue());
			Map<String, Object> attr = new LinkedHashMap<String, Object>(4);
			if (lc != null) {
				attr.put("NC", lc.getNucCount());
				if (lc.getNucCount() > 0) {
					attr.put("NE", lc.getNucEnt());
					attr.put("GC", lc.getGcBias());
				}
				attr.put("TC", lc.getTrinucCount());
				if (lc.getTrinucCount() > 0) {
					attr.put("TE", lc.getTriEnt());
					if (lc.getTrinucCount() > 1) attr.put("DU", lc.getDUST());
				}
				switch (whence) {
				case CENTER:
					attr.put("ST",lc.getLocus().getStart());
				case START:  				
					attr.put("ED", lc.getLocus().getStart() + lc.size() - 1);  break;
				case END:
					attr.put("ST",lc.getLocus().getStart());
				}
				if (lc.getLocus().compareTo(loc) < 0)
					loc = lc.getLocus();
			}
			Genotype gt = new Genotype(gws.getKey(), noCall, 1.0, noFilters,
					attr, false);
			gc.add(gt);
		}
		vcb.genotypes(gc);
		vcb.loc(loc);
		VariantContext vc = vcb.make();
		writer.add(vc);
	}

	@Override
	public boolean isReduceByInterval() {
		return true;
	}

	public MultiWindowSequenceComplexity reduce(ReferenceContext rc,
			MultiWindowSequenceComplexity sum) {
		List<Map<Integer, SequenceComplexity.LocusComplexity>> lcm = sum.count(
				rc);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : lcm)
			if (l.size() != 0)
				emit(l);
		return sum;
	}

}
