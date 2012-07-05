package net.malariagen.gatk.coverage;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;


import net.malariagen.gatk.annotators.FragmentStartCount;
import net.malariagen.gatk.annotators.UniquenessScore;
import net.sf.samtools.SAMReadGroupRecord;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.annotator.DepthOfCoverage;
import org.broadinstitute.sting.gatk.walkers.annotator.MappingQualityZero;
import org.broadinstitute.sting.gatk.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

@By(DataSource.REFERENCE)
@PartitionBy(PartitionType.CONTIG)
@ReadFilters()
public class CoverageBiasWalker extends
		LocusWalker<CoverageBiasContext, CoverageBiasCovariateCounts> implements
		AnnotatorCompatibleWalker, TreeReducible<CoverageBiasCovariateCounts> {

	static class LocusBias {
		public ReferenceContext refContext;
		public VariantContext variantContext;
		public VariantContext forwardComplexity;
		public VariantContext reverseComplexity;
	}

	public enum GroupBy {
		SM, RG, SMRG, WS, NONE;

		public boolean implies(GroupBy o) {
			if (o == null)
				throw new IllegalArgumentException(
						"other group-by cannot be null");
			if (o == this)
				return true;
			if (this == SMRG && o != NONE && o != WS)
				return true;
			return false;
		}
	}

	@Argument(shortName = "mmq0", doc = "Maximum fraction of reads with MQ = 0 for a site to be counted", required = false)
	public double maximumMapQual0Fraction = 1.0;

	@Argument(shortName = "uqn", doc = "Uniqueness file, in  provided the uniquenes score value will be added to the output", required = false)
	public RodBinding<Feature> uniqueness = null;

	@Argument(shortName = "fsmmq", fullName = "fragmentStartMinimumMappingQuality", doc = "minimum quality for a read to be considered for the fragment-start count outputs", required = false)
	public int fsmmq = 0;

	@Argument(shortName = "groupBy", fullName = "genotypeGroupBy", doc = "Indicates wheter there should be one genotype column per sample (SM), read-group (RG), both (SMRG) or none (NONE)", required = false)
	public GroupBy groupBy = GroupBy.RG;

	@Argument(shortName = "W", fullName = "windowSize", required = false, doc = "window-size to get stats on")
	protected Integer windowSize = null;

	@Argument(shortName = "fs", fullName = "fragmentLenghts", required = false, doc = "fragment-length output for the samples under analysis. A genotype column will be output with the corresponding median window size")
	protected File fragmentLengthsFile = null;

	// @Argument(shortName = "complexity", doc =
	// "VCF containing complexity stats per position on the reference (as produced by ReferenceComplexity tool)",
	// required = true)
	// public RodBinding<VariantContext> complexity;

	@Output(shortName = "vo", fullName = "vcfOutput", doc = "Name of the output vcf file", required = false)
	public VCFWriter vcfOutput;

	@Output(shortName = "o", fullName = "output", doc = "name of the directory where to generate the covariate count outputs", required = true)
	public File outDir;

	private Set<String> groupNames;


	private Set<InfoFieldAnnotation> annotations;

	private final SortedMap<GenomeLoc, Complexity> complexityBuffer = new TreeMap<GenomeLoc, Complexity>();

	private ReferenceComplexityWalker complexityWalker = new ReferenceComplexityWalker();

	@Argument(shortName = "mdp", fullName = "minimumDepthOfCoverage", doc = "minimum depth of coverage in a sample for a site to be considered in the counts", required = false)
	public int minimumDepthOfCoverage = 1;

	@Override
	public void initialize() {
		super.initialize();
		if (fsmmq < 0)
			throw new UserException(
					"the fragment start minimum mapping quality cannot be less than 0");
		groupNames = buildGroupNames();
		annotations = buildAnnotations();
		if (groupBy == GroupBy.WS)
			throw new UserException(
					"grouping by window size (-groupBy WS) is not supported in this walker");
		// TODO need to put this functionality in a common place, rather than
		// borrow from
		// a different walker:
		FragmentLengthsWalker.checkOutDir(outDir);
		// .

		if (vcfOutput != null) {
			Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
			for (InfoFieldAnnotation a : annotations)
				headerLines.addAll(a.getDescriptions());
			for (InfoFieldAnnotation a : annotations)
				if (!(a instanceof UniquenessScore))
					headerLines.addAll(convertInfoToFormatDescriptions(a
							.getDescriptions()));
			VCFHeader header = new VCFHeader(headerLines, groupNames);
			vcfOutput.writeHeader(header);
		}
		initializeComplexityWalker();
	}

	private void initializeComplexityWalker() {
		final GenomeLocParser locParser = this.getToolkit()
				.getGenomeLocParser();
		complexityWalker.setToolkit(this.getToolkit());
		complexityWalker.groupBy = this.groupBy;
		complexityWalker.rounding = 1;
		complexityWalker.windowSize = this.windowSize;
		complexityWalker.fragmentLengthsFile = this.fragmentLengthsFile;
		complexityWalker.writer = new VCFWriter() {

			@Override
			public void writeHeader(VCFHeader header) {
				// Ignore it.
			}

			@Override
			public void close() {
				// Ignore it.
			}

			@Override
			public void add(VariantContext vc) {
				
				GenomeLoc loc = locParser.createGenomeLoc(vc.getChr(),
						vc.getStart(), vc.getStart());
				Complexity c = complexityBuffer.get(loc);
				if (c == null) {
					complexityBuffer.put(loc, c = new Complexity(loc, vc));
				}
				else 
					c.forward = vc;
				for (Genotype gt : vc.getGenotypes()) {
					int end = gt.getAttributeAsInt("ED", -1);				
					if (end == -1)
						continue;
					GenomeLoc rloc = locParser.createGenomeLoc(loc.getContig(),
							end, end);
					Complexity rc = complexityBuffer.get(rloc);
					if (rc == null)
						complexityBuffer.put(rloc,
								rc = new Complexity(rloc, gt));
					else
						rc.addReverse(gt);
				}
			}

		};
		complexityWalker.initialize();

	}

	@Override
	public void onTraversalDone(CoverageBiasCovariateCounts sum) {
		super.onTraversalDone(sum);
		try {
			sum.saveIn(outDir);
		} catch (IOException e) {
			throw new StingException(
					"could not save the counters due to and error", e);
		}
	}

	private Collection<VCFFormatHeaderLine> convertInfoToFormatDescriptions(
			List<VCFInfoHeaderLine> descriptions) {
		List<VCFFormatHeaderLine> result = new LinkedList<VCFFormatHeaderLine>();
		for (VCFInfoHeaderLine original : descriptions)
			result.add(new VCFFormatHeaderLine(
					original.getID(),
					original.getCountType(),
					original.getType(),
					"Same as info field with the same ID but stratified per sample and or read-group"));
		return result;
	}

	private Set<InfoFieldAnnotation> buildAnnotations() {
		Set<InfoFieldAnnotation> result = new LinkedHashSet<InfoFieldAnnotation>();
		result.add(new RMSMappingQuality());
		result.add(new MappingQualityZero());
		if (uniqueness != null)
			result.add(new UniquenessScore());
		result.add(new DepthOfCoverage());
		result.add(new FragmentStartCount(fsmmq));
		return result;
	}

	private Set<String> buildGroupNames() {
		Set<String> result = new LinkedHashSet<String>();
		if (groupBy.implies(GroupBy.SM)) // so includes SM
			for (Sample s : getToolkit().getSampleDB().getSamples())
				if (!result.add(s.getID()))
					throw new UserException(
							"there is more than one group or sample that share id "
									+ s.getID());
		if (groupBy.implies(GroupBy.RG)) // so includes RG
			for (SAMReadGroupRecord rg : getToolkit().getReadsDataSource()
					.getHeader().getReadGroups())
				if (!result.add(rg.getId()))
					throw new UserException(
							"there is more than one group or sample that share id "
									+ rg.getId());
		return result;
	}

	@Override
	public CoverageBiasCovariateCounts reduceInit() {
		return new CoverageBiasCovariateCounts(groupNames,
				complexityWalker.reduceInit());
	}

	@Override
	public CoverageBiasCovariateCounts reduce(CoverageBiasContext value,
			CoverageBiasCovariateCounts sum) {

		sum.complexity = complexityWalker.reduce(value.getReferenceContext(),
				sum.complexity);

		if (vcfOutput != null)
			vcfOutput.add(value);

		GenomeLoc loc = value.getReferenceContext().getLocus();
		GenomeLoc toLoc = getToolkit().getGenomeLocParser().incPos(loc);

		Map<GenomeLoc, Complexity> headBuffer = complexityBuffer.headMap(toLoc);
		GenomeLoc stopLoc = null;
		for (Complexity c : headBuffer.values()) {
			VariantContext forwardComplexityVc = c.forward;
			if (c.forward == null) {
				stopLoc = c.locus;
				break;
			}
			VariantContext reverseComplexityVc = new VariantContextBuilder(
					c.forward).genotypes(c.reverse).make();
			GenotypesContext gc = value.getGenotypes();
			for (Genotype g : gc) {
				Genotype forwardComplexity = forwardComplexityVc == null ? null
						: forwardComplexityVc.getGenotype(g.getSampleName());
				Genotype reverseComplexity = reverseComplexityVc == null ? null
						: reverseComplexityVc.getGenotype(g.getSampleName());

				double fGcBias = forwardComplexity == null ? Double.NaN
						: forwardComplexity.getAttributeAsDouble("GC",
								Double.NaN);
				double rGcBias = reverseComplexity == null ? Double.NaN
						: reverseComplexity.getAttributeAsDouble("GC",
								Double.NaN);

				if (Double.isNaN(rGcBias) && Double.isNaN(fGcBias))
					continue;
				
				int size = forwardComplexity != null ? forwardComplexity.getAttributeAsInt("ED",c.locus.getStart() + 1) - c.locus.getStart() : 1;
				int dp = g.getAttributeAsInt("DP", -1);
				if (dp < minimumDepthOfCoverage)
					continue;
				if (maximumMapQual0Fraction < 1.0) {
					int mq0 = g.getAttributeAsInt("MQ0",-1);
					if (maximumMapQual0Fraction * dp < mq0)
						continue;
				}
				int ffs = g.getAttributeAsInt(FragmentStartCount.FORWARD_START_KEY,0);
				int rfs = g.getAttributeAsInt(FragmentStartCount.REVERSE_START_KEY,0);
				if (!Double.isNaN(fGcBias) && size > 0)
					sum.add(g.getSampleName(), fGcBias, size, 1, ffs,true);
				if (!Double.isNaN(rGcBias) && size > 0)
					sum.add(g.getSampleName(), rGcBias, size, 1, rfs,false);
			}
		}
		if (stopLoc == null)
			headBuffer.clear();
		else
			complexityBuffer.headMap(stopLoc).clear();
		return sum;
	}

	@Override
	public CoverageBiasCovariateCounts treeReduce(
			CoverageBiasCovariateCounts lhs, CoverageBiasCovariateCounts rhs) {
		lhs.mergeIn(rhs);
		return lhs;
	}

	@Override
	public CoverageBiasContext map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context) {
		complexityWalker.map(tracker, ref, context);
		VariantContextBuilder vcb = new VariantContextBuilder();
		vcb.loc(ref.getLocus());
		vcb.alleles(Collections.singletonList(Allele.create(ref.getBase(), true)));
		VariantContext fooVc = vcb.make();
		Map<String, Object> vcAttributes = new LinkedHashMap<String, Object>();
		Map<String, AlignmentContext> stratifiedCtx = stratifyByGroupName(context);
		for (InfoFieldAnnotation anno : annotations)
			vcAttributes.putAll(anno.annotate(tracker, this, ref,
					stratifiedCtx, fooVc));
		vcb.attributes(vcAttributes);
		if (groupBy != GroupBy.NONE)
			vcb.genotypes(buildGenotypes(tracker, ref, fooVc, stratifiedCtx));
		return new CoverageBiasContext(vcb.make(), ref);
	}

	private GenotypesContext buildGenotypes(RefMetaDataTracker tracker,
			ReferenceContext ref, VariantContext fooVc,
			Map<String, AlignmentContext> stratifiedCtx) {
		GenotypesContext gc = GenotypesContext.create();
		Set<String> noFilters = Collections.emptySet();
		List<Allele> noCalls = Collections.singletonList(Allele.NO_CALL);
		for (String groupName : groupNames) {
			Map<String, Object> gAttributes = new LinkedHashMap<String, Object>(
					10);
			Map<String, AlignmentContext> fooCtx = Collections.singletonMap(
					groupName, stratifiedCtx.get(groupName));
			for (InfoFieldAnnotation anno : annotations)
				if (!(anno instanceof UniquenessScore))
					gAttributes.putAll(anno.annotate(tracker, this, ref,
							fooCtx, fooVc));
			Genotype gt = new Genotype(groupName, noCalls, 1, noFilters,
					gAttributes, false);

			gc.add(gt);
		}
		return gc;
	}

	private Map<String, AlignmentContext> stratifyByGroupName(
			AlignmentContext context) {
		if (groupNames.size() == 0) // i.e. groupBy == NONE
			return Collections.singletonMap("ALL", context);
		Map<String, List<PileupElement>> lists = new LinkedHashMap<String, List<PileupElement>>(
				groupNames.size());
		for (String groupName : groupNames)
			lists.put(groupName, new ArrayList<PileupElement>(context.size()));

		for (PileupElement pe : context.getBasePileup()) {
			SAMReadGroupRecord rg = pe.getRead().getReadGroup();
			boolean groupFound = false;
			if (groupBy.implies(GroupBy.RG) && groupNames.contains(rg.getId())) {
				lists.get(rg.getId()).add(pe);
				groupFound = true;
			}
			if (groupBy.implies(GroupBy.SM)) {
				String s = rg.getSample();
				if (groupNames.contains(s)) {
					lists.get(s).add(pe);
					groupFound = true;
				}
			}
			if (!groupFound)
				throw new IllegalArgumentException(
						"there is some orphan read in context "
								+ pe.getRead().getReadName());

		}
		Map<String, AlignmentContext> result = new LinkedHashMap<String, AlignmentContext>(
				groupNames.size());
		for (String groupName : groupNames)
			result.put(groupName,
					new AlignmentContext(context.getLocation(),
							new ReadBackedPileupImpl(context.getLocation(),
									lists.get(groupName))));
		return result;
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

	private class Complexity {

		public final GenomeLoc locus;
		public VariantContext forward;
		public GenotypesContext reverse;

		Complexity(GenomeLoc l, VariantContext fwd) {
			locus = l;
			forward = fwd;
			reverse = GenotypesContext.create();
		}

		Complexity(GenomeLoc l, Genotype gt) {
			this(l, (VariantContext) null);
			reverse.add(gt);
		}

		void addReverse(Genotype gt) {
			reverse.add(gt);
		}

		double getGcBias(String name, boolean forward) {
			if (forward)
				return getForwardGcBias(name);
			else
				return getReverseGcBias(name);
		}

		private double getReverseGcBias(String name) {
			Genotype gt = reverse.get(name);
			if (gt == null)
				return Double.NaN;
			return gt.getAttributeAsDouble("GC", Double.NaN);
		}

		private double getForwardGcBias(String name) {
			if (forward == null)
				return Double.NaN;
			Genotype gt = forward.getGenotype(name);
			if (gt == null)
				return Double.NaN;
			return gt.getAttributeAsDouble("GC", Double.NaN);
		}
	}

}
