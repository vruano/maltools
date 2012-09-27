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
import net.malariagen.gatk.gff.GFFFeature;
import net.malariagen.gatk.uniqueness.UQNFeature;
import net.malariagen.gatk.walker.SequenceComplexity.LocusComplexity;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
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
@Requires({DataSource.REFERENCE,DataSource.REFERENCE_BASES,DataSource.REFERENCE_ORDERED_DATA})
@PartitionBy(PartitionType.CONTIG)
@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
public class ReferenceComplexityWalker
		extends
		LocusWalker<net.malariagen.gatk.walker.ReferenceComplexityWalker.Locus, MultiWindowSequenceComplexity>
		implements AnnotatorCompatibleWalker {

	public enum Whence {
		START, CENTER, END;
	}

	static class Locus {
		ReferenceContext ref;
		public int uniqueness;
		public boolean coding;
	}

	@Argument(shortName = "W", fullName = "windowSize", required = false, doc = "window-size to get stats on")
	protected List<Integer> windowSizeList = Collections.emptyList();

	@Argument(shortName = "useInfoFields", fullName = "useInfoFields", required = false, doc = "whether complexity for a fixed unique window size shoul be outputted in the info fields")
	protected boolean useInfoFields = false;

	protected Integer windowSize;

	@Argument(shortName = "uniqueness", doc = "Uniqueness file, in  provided the uniquenes score value will be added to the output", required = false)
	public RodBinding<UQNFeature> uniqueness = null;

	@Argument(shortName = "features", doc = "Annotation GFF file", required = false)
	public RodBinding<GFFFeature> features = null;

	@Argument(shortName = "rounding", fullName = "windowSizeRounding", doc = "window size rounding, 1 means no rounding", required = false)
	protected int rounding = 1;

	@Output(shortName = "o", fullName = "output", doc = "File to which variants should be written", required = true)
	protected VCFWriter writer = null;

	@Argument(shortName = "whence", fullName = "whence", doc = "Where the targeted position in the reference is located within the window (START,CENTER or END) ", required = false)
	protected Whence whence = Whence.CENTER;

	private Map<GenomeLoc, Locus> locus;

	@Override
	public Locus map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		Locus result = new Locus();
		result.ref = ref;
		if (uniqueness.isBound())
			result.uniqueness = uniquenessScore(tracker, ref);
		if (features.isBound())
			result.coding = isCoding(tracker);
		locus.put(ref.getLocus(), result);
		return result;
	}

	private boolean isCoding(RefMetaDataTracker tracker) {
		for (GFFFeature gf : tracker.getValues(features)) {
			if (gf.getType().isProteinCoding())
				return true;
		}
		return false;
	}

	private int uniquenessScore(RefMetaDataTracker tracker, ReferenceContext ref) {
		UQNFeature gf = tracker.getFirstValue(uniqueness);
		return gf == null ? 99 : gf.getScore();
	}

	@Output(shortName = "groupBy", required = false)
	protected GroupBy groupBy = null;

	Map<String, Integer> groupWindowSize;

	@Override
	public void initialize() {
		super.initialize();
		locus = new HashMap<GenomeLoc, Locus>();
		if (writer == null)
			throw new IllegalStateException("the write is yet null");
		if (!windowSizeList.isEmpty())
			windowSize = windowSizeList.get(0);
		else
			throw new IllegalArgumentException(
					"you must indicate at least one window size (-W)");
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
				throw new UserException(
						"no requested window size can be null, zero or negative (-W "
								+ ws + ")");
			int actualWs = (int) Math.round(ws / rounding) * rounding;
			groupWindowSize.put("" + ws, actualWs);
		}
		Set<String> groupNames;
		if (uniqueness.isBound())
			headerLines.add(new VCFInfoHeaderLine("UQ", 1,
					VCFHeaderLineType.Integer, "Uniqueness score"));
		if (features.isBound())
			headerLines.add(new VCFInfoHeaderLine("CODING", 0,
					VCFHeaderLineType.Flag,
					"Marks positions that are withing protein coding exons"));
		if (useInfoFields) {
			groupNames = Collections.emptySet();
			if (features.isBound())
				headerLines
						.add(new VCFInfoHeaderLine("CodingCount", 1,
								VCFHeaderLineType.Integer,
								"Number of positions that are considered coding within the window; can be greater than NC if the window contains N or X"));
			headerLines
					.add(new VCFInfoHeaderLine("GCCnt", 1,
							VCFHeaderLineType.Integer,
							"Number of G or C nucleotides in the window"));
			headerLines.add(new VCFInfoHeaderLine("NucEnt", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy in nats (only if NucCnt > 0)"));
			headerLines.add(new VCFInfoHeaderLine("TriEnt", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy in nats (only if TriCnt > 0)"));
			headerLines
					.add(new VCFInfoHeaderLine("NucCnt", 1,
							VCFHeaderLineType.Integer,
							"Number of sense nucleotides in the window (A, C, G or T)"));
			headerLines
					.add(new VCFInfoHeaderLine("TriCnt", 1,
							VCFHeaderLineType.Integer,
							"Number of sense trinucleotides in the window (composed only by A, C, G or T)"));
			headerLines
					.add(new VCFInfoHeaderLine("DUST", 1,
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

		} else {
			groupNames = groupWindowSize.keySet();
			if (features.isBound())
				headerLines
						.add(new VCFFormatHeaderLine("CC", 1,
								VCFHeaderLineType.Integer,
								"Number of positions that are considered coding within the window; can be greater than NC if the window contains N or X"));
			headerLines
					.add(new VCFFormatHeaderLine("GT", 1,
							VCFHeaderLineType.String,
							"Genotype Call, never used but left in to be compliant with all standards"));
			headerLines
					.add(new VCFFormatHeaderLine("GC", 1,
							VCFHeaderLineType.Integer,
							"Number of G or C nucleotides in the window"));
			headerLines.add(new VCFFormatHeaderLine("NE", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy in nats (only if NC > 0)"));
			headerLines.add(new VCFFormatHeaderLine("TE", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy in nats (only if TC > 0)"));
			headerLines
					.add(new VCFFormatHeaderLine("NC", 1,
							VCFHeaderLineType.Integer,
							"Number of sense nucleotides in the window (A, C, G or T)"));
			headerLines
					.add(new VCFFormatHeaderLine("TC", 1,
							VCFHeaderLineType.Integer,
							"Number of sense trinucleotides in the window (composed only by A, C, G or T)"));
			headerLines
					.add(new VCFFormatHeaderLine("DU", 1,
							VCFHeaderLineType.Float,
							"DUST complexity metric for the given window (only if TC > 1)"));
			headerLines.add(new VCFFormatHeaderLine("ED", 1,
					VCFHeaderLineType.Float,
					"Stop position of the window interval"));
			switch (whence) {
			case CENTER:
				headerLines.add(new VCFFormatHeaderLine("ST", 1,
						VCFHeaderLineType.Float,
						"Start position of the window interval"));
			case START:
				headerLines.add(new VCFFormatHeaderLine("ED", 1,
						VCFHeaderLineType.Float,
						"Stop position of the window interval"));
				break;
			case END:
				headerLines.add(new VCFFormatHeaderLine("ST", 1,
						VCFHeaderLineType.Float,
						"Start position of the window interval"));
			}

		}
		header = new VCFHeader(headerLines, groupNames);
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
		return new MultiWindowSequenceComplexity(windowSizeInts, whence, this
				.getToolkit().getGenomeLocParser());
	}

	@Override
	public MultiWindowSequenceComplexity reduce(Locus value,
			MultiWindowSequenceComplexity sum) {
		if (sum.lastLocus() != null
				&& value.ref.getLocus().compareContigs(sum.lastLocus()) != 0) {
			for (Map<Integer, SequenceComplexity.LocusComplexity> l : sum
					.flush())
				if (l.size() != 0)
					emit(l);
			sum = this.reduceInit();
		}
		List<Map<Integer, SequenceComplexity.LocusComplexity>> lcm = sum.count(
				value.ref, value.coding);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : lcm)
			if (l.size() != 0)
				emit(l);
		return sum;
	}

	@Override
	public void onTraversalDone(MultiWindowSequenceComplexity sum) {
		super.onTraversalDone(sum);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : sum.flush())
			if (l.size() != 0)
				emit(l);
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
			if (features.isBound())
				attributes.put("CodingCount", lc.getCodingCount());
			if (lc.getNucCount() > 0) {
				attributes.put("NucEnt", lc.getNucEnt());
			}
			if (lc.getTrinucCount() > 0) {
				attributes.put("TriEnt", lc.getTriEnt());
				if (lc.getTrinucCount() > 1)
					attributes.put("DUST", lc.getDUST());
			}
			attributes.put("GCCnt", lc.getGcCount());
			attributes.put("NucCnt", lc.getNucCount());
			attributes.put("TriCnt", lc.getTrinucCount());
			switch (whence) {
			case CENTER:
				attributes.put("START", lc.getLocus().getStart());
			case START:
				attributes.put("END", lc.getLocus().getStart() + lc.size() - 1);
				break;
			case END:
				attributes.put("START", lc.getLocus().getStart());
			}
			if (lc.getLocus().compareTo(loc) < 0)
				loc = lc.getLocus();
		}
		Locus l = locus.remove(loc);
		if (uniqueness.isBound())
			attributes.put("UQ", l.uniqueness);
		if (features.isBound() && l.coding)
			attributes.put("CODING", null);
		vcb.attributes(attributes);
		List<Allele> noCall = Collections.singletonList(Allele.NO_CALL);
		Set<String> noFilters = Collections.emptySet();
		GenotypesContext gc = GenotypesContext.create();
		for (Map.Entry<String, Integer> gws : groupWindowSize.entrySet()) {
			LocusComplexity lc = lcm.get(gws.getValue());
			Map<String, Object> attr = new LinkedHashMap<String, Object>(4);
			if (lc != null) {
				if (features.isBound())
					attr.put("CC", lc.getCodingCount());
				attr.put("NC", lc.getNucCount());
				attr.put("GC", lc.getGcCount());
				if (lc.getNucCount() > 0) {
					attr.put("NE", lc.getNucEnt());
				}
				attr.put("TC", lc.getTrinucCount());
				if (lc.getTrinucCount() > 0) {
					attr.put("TE", lc.getTriEnt());
					if (lc.getTrinucCount() > 1)
						attr.put("DU", lc.getDUST());
				}
				switch (whence) {
				case CENTER:
					attr.put("ST", lc.getLocus().getStart());
				case START:
					attr.put("ED", lc.getLocus().getStart() + lc.size() - 1);
					break;
				case END:
					attr.put("ST", lc.getLocus().getStart());
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
		List<Map<Integer, SequenceComplexity.LocusComplexity>> lcm = sum
				.count(rc);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : lcm)
			if (l.size() != 0)
				emit(l);
		return sum;
	}

	@Override
	public RodBinding<VariantContext> getSnpEffRodBinding() {
		return null;
	}

	@Override
	public RodBinding<VariantContext> getDbsnpRodBinding() {
		return null;
	}

	@Override
	public List<RodBinding<VariantContext>> getCompRodBindings() {
		return Collections.emptyList();
	}

	@Override
	public List<RodBinding<VariantContext>> getResourceRodBindings() {
		return Collections.emptyList();
	}

	@Override
	public boolean alwaysAppendDbsnpId() {
		return false;
	}

}
