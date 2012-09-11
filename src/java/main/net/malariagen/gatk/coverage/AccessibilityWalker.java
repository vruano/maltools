package net.malariagen.gatk.coverage;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.coverage.CoverageBiasWalker.GroupBy;
import net.malariagen.gatk.gff.GFFFeature;
import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.math.IntegerDistributionSet;
import net.malariagen.gatk.utils.ReadGroupDB;
import net.malariagen.utils.Nucleotide;
import net.sf.picard.reference.IndexedFastaSequenceFile;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode;
import org.broadinstitute.sting.utils.baq.BAQ.QualityMode;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElementFilter;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

@ReadFilters({NotPrimaryAlignmentFilter.class,UnmappedReadFilter.class})
@PartitionBy(PartitionType.LOCUS)
@By(DataSource.REFERENCE)
@Requires({DataSource.REFERENCE_BASES, DataSource.READS, DataSource.REFERENCE_ORDERED_DATA})
public class AccessibilityWalker extends LocusWalker<VariantContext, AccessibilityStatistics> implements
TreeReducible<AccessibilityStatistics>{

	@Argument(fullName = "minimumBaseQuality", shortName = "mbq", doc = "Minimum quality for a base to be considered in coverage counting, -1 to indicate no minimum (default)", required = false)
	public int minimumBaseQuality = -1;

	@Argument(fullName = "minimumMappingQuality", shortName = "mmq", doc = "Minimum mapping quality for a read (thus its based) to be considered in coverage counting, -1 to indicate no minimum (default)", required = false)
	public int minimumMappingQuality = -1;

	@Argument(fullName = "minimumBAQ", shortName = "mbaq", doc = "Minimum BAQ for a base to be considered in coverage counting, -1 to indicate no minimum (default)", required = false)
	public int minimumBAQ = -1;

	@Argument(fullName = "excludeAmbigousRef", shortName = "exclAmbRef", doc = "Indicates whether totally degenerated ambiguous reference sites should be excluded", required = true)
	public boolean excludeAmbigousRef = false;

	@Argument(fullName = "excludeAmbigousCall", shortName = "exclAmbCall", doc = "Indicates whether totally degenerated ambiguous read-bases should be excluded", required = false)
	public boolean excludeAmbigousBase = false;

	@Argument(fullName = "includeReadDeletions", shortName = "inclRdDel", doc = "Do not consider read deletions in coverage depth calculation", required = false)
	public boolean includeReadDeletions = false;
		
	@Argument(fullName="normalizePerChromosome", doc="whether we should normalize for each chromosome separately (implies -normalize)", required=false)
	protected boolean normalizePerChromosome = false;
	
	@Argument(fullName="normalizePerCodingStatus", doc="whether we should normalize per coding-status (implies -normalize)", required=false)
	protected boolean normalizePerCodingStatus = false;
	
	@Argument(fullName="normalize", shortName="normalize", doc="whether we should normalize by median depth", required=false)
	protected Boolean normalize = null;
	
	@Argument(fullName="features", shortName="features", doc="file containing the GFF records indicating where the coding regions are located", required=false)
	protected RodBinding<GFFFeature> features = null;
	
	@Argument(fullName="groupBy", shortName="groupBy", doc="whether we should generate stats per read-group (RG), sample (SM), both (SMRG) or none (NONE default)", required=false)
	protected GroupBy groupBy = GroupBy.NONE;
	
	protected ReadGroupDB readGroupDb;
	
	@Argument(shortName="cds", fullName="coverageDistributionSet", doc="coverage distribution set (implies -normalize)", required=false)
	protected List<File> coverageDistributionFiles = Collections.emptyList();

	@Output(shortName="o", doc="output vcf file name, by default it uses the standard output", required=false)
	protected VCFWriter writer;
	
	IntegerDistributionSet coverageDistributionSet;

    Set<String> groupNames;
    
    private ThreadLocal<AccesibilityAnnotations> annotations;

	private PileupElementFilter pileupFilter;
	
	@Override
	public void initialize() {
		super.initialize();
		
		if (!coverageDistributionFiles.isEmpty()) {
			List<IntegerDistributionSet> cdsList = new LinkedList<IntegerDistributionSet>();
			for (File f : coverageDistributionFiles) {
				if (f == null)
					throw new IllegalArgumentException("coverage-distribution-files cannot contain nulls");
				if (!f.isFile())
					throw new UserException("coverage-distribution file (" + f + ") is unreachable or is not a regular file");
				if (!f.canRead())
					throw new UserException("coverage-distribution file (" + f + ") cannot be read");
				try {
				  Reader r = new FileReader(f);
				  cdsList.add(IntegerDistributionSet.read(r));
				}
				catch (IOException e) {
				  throw new StingException("unexpected exception");	
				}
			}
			coverageDistributionSet = IntegerDistributionSet.merge(cdsList.toArray(new IntegerDistributionSet[cdsList.size()]));
		}
		
		readGroupDb = new ReadGroupDB(getToolkit());
		groupNames = new HashSet<String>();
		
		if (groupBy == GroupBy.WS)
			throw new UserException("group-by WS is not applicable in this walker");
		if (groupBy.implies(GroupBy.SM))
			groupNames.addAll(readGroupDb.getSampleIDs());
		if (groupBy.implies(GroupBy.RG))
			groupNames.addAll(readGroupDb.getReadGroupIDs());
		
		
		if (normalize == null && (normalizePerCodingStatus || normalizePerChromosome))
			normalize = true;
		else if (normalize != null && normalize == false) {
			if (normalizePerCodingStatus || normalizePerChromosome)
				throw new UserException("incoherent normalization flag argument set-up");
		}
		else if (normalize != null && normalize) {
			if (coverageDistributionSet == null)
				throw new UserException("you cannot request normalization without providing at least one coverage distribution set file");
		}
		else
			normalize = coverageDistributionSet != null;
		
		if (normalize) {
			for (String gn : groupNames)
				if (coverageDistributionSet.getSampleDistributionSet(gn) == null)
					throw new UserException("the coverage distribution set provided does not contain information for group '" + gn + "'");
			if (normalizePerChromosome) 
				for (GenomeLoc loc : getToolkit().getIntervals())
				if (coverageDistributionSet.getAllSamplesDistributionSet().getSequenceDistributionSet(loc.getContig()) == null)
					throw new UserException("the coverage distribution set provided does not contain information for chromosome/sequence '" + loc.getContig() + "'");
			if (normalizePerCodingStatus) {
				if (features == null)
					throw new UserException("per coding-status normalization requested but no genome feature GFF file provided (-features option) ");
				IntegerDistribution codingDist = coverageDistributionSet.getAllSamplesDistributionSet().getSequenceDistributionSet(getToolkit().getIntervals().iterator().next().getContig()).getCategoryDistribution(LocusCategory.CODING);
				IntegerDistribution nonCodingDist = coverageDistributionSet.getAllSamplesDistributionSet().getSequenceDistributionSet(getToolkit().getIntervals().iterator().next().getContig()).getCategoryDistribution(LocusCategory.NON_CODING);
				
				if (codingDist == null || codingDist.count() == 0)
					throw new UserException("coding status normalization request but the coverage-distribution set provided does not seem to contain such information");
				if (nonCodingDist == null || nonCodingDist.count() == 0)
					throw new UserException("coding status normalization request but the coverage-distribution set provided does not seem to contain such information");
			}
		}
		
		annotations = new ThreadLocal<AccesibilityAnnotations>();
		initializeBAQCalculatingEngine();
		pileupFilter = CoverageDistributionWalker.initializePileupFilter(!includeReadDeletions, excludeAmbigousBase, minimumBaseQuality, minimumMappingQuality, minimumBAQ, baqHMM, reference, cmode, qmode);

		writer.writeHeader(new VCFHeader(AccesibilityAnnotations.getHeaderLines(this),groupNames));
	}
	
	private BAQ baqHMM;
	private IndexedFastaSequenceFile reference;
	private CalculationMode cmode;
	private QualityMode qmode;

	private void initializeBAQCalculatingEngine() {
		baqHMM = new BAQ();
		reference = getToolkit().getReferenceDataSource().getReference();
		qmode = getToolkit().getWalkerBAQQualityMode();
		cmode = getToolkit().getArguments().BAQMode;
	}

	@Override
	public VariantContext map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context) {
		if (excludeAmbigousRef && !BaseUtils.isRegularBase(ref.getBase()))
			return null;
		Nucleotide refNuc = Nucleotide.fromByte(ref.getBase());
		List<Allele> alleles = new LinkedList<Allele>();
		Map<String, AlignmentContext> stratified = CoverageBiasWalker.stratifyByGroupName(context, groupNames, groupBy);
		ReadBackedPileup pileup = context.getBasePileup().getFilteredPileup(
				pileupFilter);
//		if (pileup.isEmpty()) 
//			return null;
		VariantContextBuilder result = new VariantContextBuilder();
		result.loc(ref.getLocus());
		result.referenceBaseForIndel(Nucleotide.N.byteValue());
		int[] baseCounts = pileup.getBaseCounts();
		Allele refAllele = Allele.create(ref.getBase(),true);
		result.alleles(alleles);
		Map<String, Object> infoAnnotations = new LinkedHashMap<String,Object>();
		Map<String, Map<String, Object>> formatAnnotations = new LinkedHashMap<String,Map<String,Object>>(stratified.size());
		for (String s : stratified.keySet())
			formatAnnotations.put(s, new LinkedHashMap<String,Object>(10));
		AccesibilityAnnotations anno = annotations.get();
		if (anno == null)
			annotations.set(anno = new AccesibilityAnnotations(this));
		
		anno.annotate(tracker, ref, stratified, infoAnnotations, formatAnnotations);
		GenotypesContext gc = GenotypesContext.create();
		result.attributes(infoAnnotations);
		Set<String> noFilters = Collections.emptySet();
		List<Allele> noCall = Collections.emptyList();
		for (String s : stratified.keySet()) {
			Genotype g = new Genotype(s,noCall,Genotype.NO_LOG10_PERROR,noFilters,formatAnnotations.get(s),false);
			gc.add(g);
		}
		result.genotypes(gc);
		alleles.add(refAllele);
		for (int i = 0; i < 4; i++)
			if (refNuc.ordinal() == i)
				continue;
			else if (baseCounts[i] != 0)
				alleles.add(Allele.create(Nucleotide.values()[i].byteValue(),false));
		return result.make();
	}


	@Override
	public AccessibilityStatistics reduceInit() {
		return new AccessibilityStatistics();
	}

	@Override
	public AccessibilityStatistics reduce(VariantContext value,
			AccessibilityStatistics sum) {
		if (value != null) writer.add(value);
		return sum;
	}

	public boolean isCoding(RefMetaDataTracker tracker) {
		for (GFFFeature gf : tracker.getValues(GFFFeature.class))
				if (gf.getType().isProteinCoding())
					return true;
		return false;
	}


	@Override
	public AccessibilityStatistics treeReduce(AccessibilityStatistics lhs,
			AccessibilityStatistics rhs) {
		return lhs;
	}

}
