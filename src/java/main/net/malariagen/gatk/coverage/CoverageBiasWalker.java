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
import java.util.concurrent.ConcurrentHashMap;

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
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.annotator.DepthOfCoverage;
import org.broadinstitute.sting.gatk.walkers.annotator.MappingQualityZero;
import org.broadinstitute.sting.gatk.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.GenomeLoc;
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
@PartitionBy(PartitionType.LOCUS)
public class CoverageBiasWalker
		extends
		LocusWalker<net.malariagen.gatk.coverage.CoverageBiasWalker.LocusBias, CoverageBiasCovariateCounts>
		implements AnnotatorCompatibleWalker,
		TreeReducible<CoverageBiasCovariateCounts> {

	static class LocusBias {
		public VariantContext variantContext;
		public VariantContext forwardComplexity;
		public VariantContext reverseComplexity;
	}

	public enum GroupBy {
		SM, RG, SMRG, NONE;

		public boolean implies(GroupBy o) {
			if (o == null)
				throw new IllegalArgumentException(
						"other group-by cannot be null");
			if (o == this)
				return true;
			if (this == SMRG && o != NONE)
				return true;
			return false;
		}
	}

	@Argument(shortName = "uqn", doc = "Uniqueness file, in  provided the uniquenes score value will be added to the output", required = false)
	public RodBinding<Feature> uniqueness = null;

	@Argument(shortName = "fsmmq", fullName = "fragmentStartMinimumMappingQuality", doc = "minimum quality for a read to be considered for the fragment-start count outputs", required = false)
	public int fsmmq = 0;

	@Argument(shortName = "groupBy", fullName = "genotypeGroupBy", doc = "Indicates wheter there should be one genotype column per sample (SM), read-group (RG), both (SMRG) or none (NONE)", required = false)
	public GroupBy groupBy = GroupBy.RG;

	@Argument(shortName = "complexity", doc = "VCF containing complexity stats per position on the reference (as produced by ReferenceComplexity tool)", required = true)
	public RodBinding<VariantContext> complexity;

	@Output(shortName = "vo", fullName = "vcfOutput", doc = "Name of the output vcf file", required = false)
	public VCFWriter vcfOutput;

	@Output(shortName = "o", fullName = "output", doc = "name of the directory where to generate the covariate count outputs", required = true)
	public File outDir;
	
	private Set<String> groupNames;

	private Map<GenomeLoc, VariantContext> reverseComplexity;

	private Set<InfoFieldAnnotation> annotations;

	@Override
	public void initialize() {
		super.initialize();
		if (fsmmq < 0)
			throw new UserException(
					"the fragment start minimum mapping quality cannot be less than 0");
		groupNames = buildGroupNames();
		annotations = buildAnnotations();
		reverseComplexity = new ConcurrentHashMap<GenomeLoc, VariantContext>(
				10000);
		//TODO need to put this functionality in a common place, rather than borrow from
		// a different walker:
		FragmentLengthWalker.checkOutDir(outDir);
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
	}
	
	@Override
	public void onTraversalDone(CoverageBiasCovariateCounts sum) {
		super.onTraversalDone(sum);
		try {
			sum.saveIn(outDir);
		} catch (IOException e) {
			throw new StingException("could not save the counters due to and error",e);
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
		return new CoverageBiasCovariateCounts(groupNames);
	}

	@Override
	public CoverageBiasCovariateCounts reduce(LocusBias value,
			CoverageBiasCovariateCounts sum) {
		if (vcfOutput != null) vcfOutput.add(value.variantContext);
		double fGcBias, fNucEnt, fTrinucEnt, rGcBias, rNucEnt, rTrinucEnt;
		if (value.forwardComplexity != null) {
			fGcBias = value.forwardComplexity.getAttributeAsDouble("GCBias",
					Double.NaN);
			fNucEnt = value.forwardComplexity.getAttributeAsDouble("NucEnt",
					Double.NaN);
			fTrinucEnt = value.forwardComplexity.getAttributeAsDouble("TriEnt",
					Double.NaN);
		} else
			fGcBias = fNucEnt = fTrinucEnt = Double.NaN;
		if (value.reverseComplexity != null) {
			rGcBias = value.reverseComplexity.getAttributeAsDouble("GCBias",
					Double.NaN);
			rNucEnt = value.reverseComplexity.getAttributeAsDouble("NucEnt",
					Double.NaN);
			rTrinucEnt = value.reverseComplexity.getAttributeAsDouble("TriEnt",
					Double.NaN);
		} else
			rGcBias = rNucEnt = rTrinucEnt = Double.NaN;
		if (Double.isNaN(rGcBias) && Double.isNaN(fGcBias))
			return sum;

		GenotypesContext gc = value.variantContext.getGenotypes();
		for (Genotype g : gc) {
			Map<String, Object> gAttr = g.getAttributes();
			Object ofs = gAttr.get(FragmentStartCount.FORWARD_START_KEY);
			int ffs = ofs instanceof Number ? ((Number) ofs).intValue()
					: Integer.parseInt(ofs.toString());
			Object ors = gAttr.get(FragmentStartCount.REVERSE_START_KEY);
			int rfs = ofs instanceof Number ? ((Number) ors).intValue()
					: Integer.parseInt(ors.toString());
			if (!Double.isNaN(fGcBias))
				sum.add(g.getSampleName(), fGcBias, fNucEnt, fTrinucEnt, 1, ffs);
			if (!Double.isNaN(rGcBias))
				sum.add(g.getSampleName(), rGcBias, rNucEnt, rTrinucEnt, 1, rfs);
		}
		return sum;
	}

	@Override
	public CoverageBiasCovariateCounts treeReduce(
			CoverageBiasCovariateCounts lhs, CoverageBiasCovariateCounts rhs) {
		lhs.mergeIn(rhs);
		return lhs;
	}

	@Override
	public LocusBias map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		VariantContext forwardComplexity = tracker
				.getFirstValue(this.complexity);
		GenomeLoc loc = ref.getLocus();
		int reversePos = forwardComplexity == null ? -1 : Integer.parseInt(forwardComplexity
				.getAttributeAsString("END", "-1"));
		if (reversePos != -1) {
			GenomeLoc revLoc = ref.getGenomeLocParser().createGenomeLoc(
					loc.getContig(), reversePos, reversePos);
			reverseComplexity.put(revLoc, forwardComplexity);
		}
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
		LocusBias result = new LocusBias();

		result.variantContext = vcb.make();
		result.forwardComplexity = forwardComplexity;
		result.reverseComplexity = reverseComplexity.remove(loc);
		return result;
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

}
