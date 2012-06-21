package net.malariagen.gatk.coverage;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.annotator.DepthOfCoverage;
import org.broadinstitute.sting.gatk.walkers.annotator.MappingQualityZero;
import org.broadinstitute.sting.gatk.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

public class MappingAnalysisWalker extends LocusWalker<VariantContext, net.malariagen.gatk.coverage.MappingAnalysisWalker.MappingStatistics> implements AnnotatorCompatibleWalker, TreeReducible<net.malariagen.gatk.coverage.MappingAnalysisWalker.MappingStatistics> {
	


	public enum GroupBy {
		SM,RG,BOTH
	}
	@Argument(shortName="U", doc="Uniqueness file, in  provided the uniquenes score value will be added to the output", required = false)
	public RodBinding<Feature> uniqueness = null;
	
	@Argument(shortName="gb", doc="Output genotype columns one per read-group", required = false)
	public GroupBy groupBy = GroupBy.SM;
	
	@Output(shortName="o", doc="Name of the output vcf file", required = true)
	public VCFWriter output;

	private Set<String> groupNames;
	
	private Set<InfoFieldAnnotation> annotations;
	
	@Override
	public void initialize() {
		super.initialize();
		Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		groupNames = buildGroupNames();
		annotations = buildAnnotations();
		for (InfoFieldAnnotation a : annotations) 
			headerLines.addAll(a.getDescriptions());
		for (InfoFieldAnnotation a : annotations)
			if (!(a instanceof UniquenessScore)) 
				headerLines.addAll(convertInfoToFormatDescriptions(a.getDescriptions()));
		VCFHeader header = new VCFHeader(headerLines,groupNames);
		output.writeHeader(header);
	}

	private Collection<VCFFormatHeaderLine> convertInfoToFormatDescriptions(
			List<VCFInfoHeaderLine> descriptions) {
		List<VCFFormatHeaderLine> result = new LinkedList<VCFFormatHeaderLine>();
		for (VCFInfoHeaderLine original : descriptions) 
			result.add(new VCFFormatHeaderLine(original.getID(), original.getCountType(), original.getType(), "Same as info field with the same ID but stratified per sample and or read-group"));
		return result;
	}

	private Set<InfoFieldAnnotation> buildAnnotations() {
		Set<InfoFieldAnnotation> result = new LinkedHashSet<InfoFieldAnnotation>();
		result.add(new RMSMappingQuality());
		result.add(new MappingQualityZero());
		result.add(new UniquenessScore());
		result.add(new DepthOfCoverage());
		result.add(new FragmentStartCount());
		return result;
	}

	private Set<String> buildGroupNames() {
		Set<String> result = new LinkedHashSet<String>();
		if (groupBy != GroupBy.RG)  // so includes SM
			for (Sample s : getToolkit().getSampleDB().getSamples())
				if (!result.add(s.getID()))
					throw new UserException("there is more than one group or sample that share id " + s.getID());
		if (groupBy != GroupBy.SM)  // so includes RG
			for (SAMReadGroupRecord rg : getToolkit().getReadsDataSource().getHeader().getReadGroups())
				if (!result.add(rg.getId()))
					throw new UserException("there is more than one group or sample that share id " + rg.getId());
		return result;
	}

	@Override
	public MappingStatistics reduceInit() {
		return new MappingStatistics();
	}

	@Override
	public MappingStatistics reduce(VariantContext value, MappingStatistics sum) {
		output.add(value);
		return sum;
	}

	@Override
	public MappingStatistics treeReduce(MappingStatistics lhs,
			MappingStatistics rhs) {
		return lhs;
	}
	
	@Override
	public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		VariantContextBuilder vcb = new VariantContextBuilder();
		vcb.alleles(Collections.singletonList(Allele.create(ref.getBase(),true)));
		VariantContext fooVc = vcb.make();
		Map<String,Object> vcAttributes = new LinkedHashMap<String,Object>();
		Map<String,AlignmentContext> stratifiedCtx = stratifyByGroupName(context);
		for (InfoFieldAnnotation anno : annotations)
			vcAttributes.putAll(anno.annotate(tracker, this, ref, stratifiedCtx, fooVc));
		
		Set<String> noFilters = Collections.emptySet();
		List<Allele> noCalls = Collections.singletonList(Allele.NO_CALL);
		GenotypesContext gc = GenotypesContext.create();
		for (String groupName : groupNames) {
			Map<String,Object> gAttributes = new LinkedHashMap<String,Object>(10);
			Map<String,AlignmentContext> fooCtx = Collections.singletonMap(groupName, stratifiedCtx.get(groupName));
			for (InfoFieldAnnotation anno : annotations)
				if (!(anno instanceof UniquenessScore))
					gAttributes.putAll(anno.annotate(tracker, this, ref, fooCtx, fooVc));
			Genotype gt = new Genotype(groupName, noCalls , 
					1, noFilters, gAttributes, false);
			gc.add(gt);
		}
		vcb.attributes(vcAttributes);
		vcb.genotypes(gc);
		return vcb.make();
	}

	private Map<String,AlignmentContext> stratifyByGroupName(
			AlignmentContext context) {
		Map<String,List<PileupElement>> lists = new LinkedHashMap<String,List<PileupElement>>(groupNames.size());
		for (String groupName : groupNames)
			lists.put(groupName, new ArrayList<PileupElement>(context.size()));
		for (PileupElement pe : context.getBasePileup()) {
			SAMReadGroupRecord rg = pe.getRead().getReadGroup();
			if (groupBy != GroupBy.SM && groupNames.contains(rg.getId())) {
				lists.get(rg.getId()).add(pe);
				continue;
			}
			if (groupBy != GroupBy.RG) {
				String s = rg.getSample();
				if (groupNames.contains(s)) {
					lists.get(s).add(pe);
					continue;
				}
			}
			throw new IllegalArgumentException("there is some orphan read in context " + pe.getRead().getReadName());
			
		}
		Map<String,AlignmentContext> result = new LinkedHashMap<String,AlignmentContext>(groupNames.size());
		for ( String groupName : groupNames ) 
			result.put(groupName, new AlignmentContext(context.getLocation(), new ReadBackedPileupImpl(context.getLocation(),lists.get(groupName))));
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
	
	public static class MappingStatistics {

	}

}
