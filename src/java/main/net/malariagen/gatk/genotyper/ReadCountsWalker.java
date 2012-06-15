package net.malariagen.gatk.genotyper;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import net.malariagen.gatk.annotators.Constants;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.variantutils.CombineVariants;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@ReadFilters({ BadMateFilter.class, UnmappedReadFilter.class })
@By(DataSource.READS)
@Downsample(by = DownsampleType.BY_SAMPLE, toCoverage = 250)
public class ReadCountsWalker extends
		LocusWalker<VariantContext, ReadCountStatistics> implements TreeReducible<ReadCountStatistics>, AnnotatorCompatibleWalker { 

	
	@Input(fullName = "coverageDistribution", shortName = "CvgD", doc = "File where to find the precalcuated coverage distribution", required = false)
	public File coverageDistributionFile;

	
	public enum UntieMethod {
		NO_CALL,
        QUALITY,
        QUALITY_AND_ORDER,
        ORDER,
	}

	public enum GenotypingMode {
		ONE_PER_ALLELE,
		MAJOR_ALLELE,
		REF_VS_ALT,
		REF_VS_ALT_01,
	}

	@Argument(fullName = "minBaseQ", shortName = "Q", doc = "minimum base quality, bases with lower quality will be discarded, 0 or lower means no filtering", required = false)
	private int minBaseQ = -1;

	@Argument(fullName = "annotation", shortName = "A", doc = "One or more specific annotations to apply to variant calls", required = false)
	protected List<String> annotationsToUse = new ArrayList<String>();

	@Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", required = false)
	protected String[] annotationClassesToUse = { "Standard" };

	@Output(fullName = "out", shortName = "o", doc = "output read counts in VCF format")
	public VCFWriter out;

	private ThreadLocal<int[]> baseToIndex = new ThreadLocal<int[]>();
	private ThreadLocal<int[]> baseCounts = new ThreadLocal<int[]>();
	private ThreadLocal<int[]> qualSums = new ThreadLocal<int[]>();
	private ThreadLocal<StringBuffer> baseCountsBuffer = new ThreadLocal<StringBuffer>();

	@Argument(fullName = "genotypingMode", shortName="gMode", required=false, doc="(def. MAJOR_ALLELE) Indicate how to produce genotypes (GT format). " + 
			                  "MAJOR_ALLELE haploid call that takes the most frequenctly seen allele, + " +
			                  "ONE_PER_ALLELE polyploid call where there is one haplo per each allele seen, " + 
			                  "REF_VS_ALT haploid call reference vs the total sume of all alternatives, major alternative is reported if the latter wins," +
			                  "REF_VS_ALT_01 haploid call same as above but always 1 (the first alterntive) is reported ending in convenient binary representation")                     
	protected GenotypingMode genotypingMode = GenotypingMode.MAJOR_ALLELE;

	
	@Argument(fullName = "untieMethod", shortName="uMethod", required=false, doc="(def. Best base quality sum) in haploid calls, how to untie when two or more reads have equal number of calls, a quality tie would result in a no-call")
	protected UntieMethod untieMethod = UntieMethod.QUALITY;	
	
	private static Allele[] refBaseToAllele = new Allele[Byte.MAX_VALUE
			- Byte.MIN_VALUE + 1];

	private static Allele[] nonRefBaseToAllele = new Allele[Byte.MAX_VALUE
			- Byte.MIN_VALUE + 1];

	private static VariantAnnotatorEngine engine = null;

	static {

		refBaseToAllele[((byte) '-') - Byte.MIN_VALUE] = refBaseToAllele[((byte) '-')
				- Byte.MIN_VALUE] = Allele.create((byte) 'A', true);
		refBaseToAllele[((byte) 'A') - Byte.MIN_VALUE] = refBaseToAllele[((byte) 'a')
				- Byte.MIN_VALUE] = Allele.create((byte) 'A', true);
		refBaseToAllele[((byte) 'C') - Byte.MIN_VALUE] = refBaseToAllele[((byte) 'c')
				- Byte.MIN_VALUE] = Allele.create((byte) 'C', true);
		refBaseToAllele[((byte) 'G') - Byte.MIN_VALUE] = refBaseToAllele[((byte) 'g')
				- Byte.MIN_VALUE] = Allele.create((byte) 'G', true);
		refBaseToAllele[((byte) 'T') - Byte.MIN_VALUE] = refBaseToAllele[((byte) 't')
				- Byte.MIN_VALUE] = Allele.create((byte) 'T', true);
		nonRefBaseToAllele[((byte) '-') - Byte.MIN_VALUE] = nonRefBaseToAllele[((byte) '-')
				- Byte.MIN_VALUE] = Allele.create((byte) '-', false);
		nonRefBaseToAllele[((byte) 'A') - Byte.MIN_VALUE] = nonRefBaseToAllele[((byte) 'a')
				- Byte.MIN_VALUE] = Allele.create((byte) 'A', false);
		nonRefBaseToAllele[((byte) 'C') - Byte.MIN_VALUE] = nonRefBaseToAllele[((byte) 'c')
				- Byte.MIN_VALUE] = Allele.create((byte) 'C', false);
		nonRefBaseToAllele[((byte) 'G') - Byte.MIN_VALUE] = nonRefBaseToAllele[((byte) 'g')
				- Byte.MIN_VALUE] = Allele.create((byte) 'G', false);
		nonRefBaseToAllele[((byte) 'T') - Byte.MIN_VALUE] = nonRefBaseToAllele[((byte) 't')
				- Byte.MIN_VALUE] = Allele.create((byte) 'T', false);
	}
	
	static {
		
		CombineVariants c;
		
		
	}

	private CoverageAnnotationEngine coverageEngine;
	
	
	@SuppressWarnings("unchecked")
	@Override
	public void initialize() {
		coverageEngine = new CoverageAnnotationEngine(annotationsToUse,coverageDistributionFile);
		coverageEngine.initializeCoverageDistribution();
		super.initialize();
		engine = new VariantAnnotatorEngine(
				Arrays.asList(annotationClassesToUse), annotationsToUse, (List<String>) Collections.EMPTY_LIST, this,getToolkit());
		Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit()
				.getSAMFileHeader());
		coverageEngine.setupCoverageAnnotation();
		out.writeHeader(new VCFHeader(getHeaderInfo(), samples));
	}

	private int[] getQualSumsArray() {
		int[] result = qualSums.get();
		if (result == null)
			qualSums.set(result = new int[20]); // just 4 possibilities ACGT but some extras for N X...
		return result;
	}	
	
	private int[] getBaseCountsArray() {
		int[] result = baseCounts.get();
		if (result == null)
			baseCounts.set(result = new int[20]); // just 4 possibilities ACGT but some extras for N X...
		return result;
	}

	private StringBuffer getBaseCountsBuffer() {
		StringBuffer result = baseCountsBuffer.get();
		if (result == null)
			baseCountsBuffer.set(result = new StringBuffer(50));
		return result;
	}

	private int[] getBaseToIndexArray() {
		int[] result = baseToIndex.get();
		if (result == null)
			baseToIndex.set(result = new int[Byte.MAX_VALUE - Byte.MIN_VALUE
					+ 1]);
		else
			Arrays.fill(result, -1);
		return result;
	}

	@Override
	public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {

		int[] baseToIndex = getBaseToIndexArray();
		int[] baseCounts = getBaseCountsArray();
		int[] qualSums = getQualSumsArray();
		StringBuffer baseCountsBuffer = getBaseCountsBuffer();

		Set<Allele> alleles = new TreeSet<Allele>();
		byte refBase = ref.getBase();
		if (refBase <= 'z' && refBase >= 'a')
			refBase = (byte) (refBase - 'a' + 'A');
		Allele refAllele = refBaseToAllele[refBase - Byte.MIN_VALUE];
		if (refAllele == null)
			refAllele = refBaseToAllele[refBase - Byte.MIN_VALUE] = Allele
					.create(refBase, true);

		alleles.add(refAllele);
		ReadBackedPileup pileup = context.getBasePileup();
		if (minBaseQ > 0)
			pileup = pileup.getBaseFilteredPileup(minBaseQ);

		for (PileupElement e : pileup) {
			byte altBase = e.getBase();
			if (altBase <= 'z' && altBase >= 'a')
				altBase = (byte) (altBase - 'a' + 'A');
			if (altBase == refBase)
				continue;
			Allele altAllele = nonRefBaseToAllele[altBase - Byte.MIN_VALUE];
			if (altAllele == null)
				altAllele = nonRefBaseToAllele[altBase - Byte.MIN_VALUE] = Allele
						.create(altBase, false);
			if (altAllele.isNull() || altAllele.isSymbolic()
					|| altAllele.isNoCall())
				continue;
			alleles.add(altAllele);
		}
		
		if (alleles.size() <= 1)
			return null;

		List<Allele> allelesList = new ArrayList<Allele>(alleles.size());
		allelesList.addAll(alleles);
		int nextIndex = 0;
		for (Allele a : allelesList) {
			byte b = a.getBases()[0];
			baseToIndex[b - Byte.MIN_VALUE] = nextIndex++;
		}

		Map<String, AlignmentContext> stratifiedContexts = AlignmentContextUtils
				.splitContextBySampleName(context);
		GenotypesContext genotypes = GenotypesContext.create();
		
		for (Map.Entry<String, AlignmentContext> e : stratifiedContexts
				.entrySet()) {
			genotypes.add(calculateGenotype(e.getKey(), e.getValue(), allelesList,
							baseToIndex, baseCounts, qualSums, baseCountsBuffer));
		}

		GenomeLoc loc = ref.getLocus();

		VariantContext vc;
		try {
			VariantContextBuilder vcb = new VariantContextBuilder();
			vcb.source("RC_call").loc(loc).alleles(allelesList).genotypes(genotypes).log10PError(1);
			vc = vcb.make();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		vc = engine.annotateContext(tracker, ref, stratifiedContexts, vc);

		return vc;
	}

	private Set<VCFHeaderLine> getHeaderInfo() {
		Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();
		headerInfo.addAll(engine.getVCFAnnotationDescriptions());

		// FORMAT and INFO fields
		headerInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
		headerInfo.add(new VCFFormatHeaderLine(Constants.READ_COUNTS_KEY,
				-1, VCFHeaderLineType.Integer,
				"Absolute read counts for each allele"));

		return headerInfo;
	}

	private Genotype calculateGenotype(String sampleName, AlignmentContext ac,
			List<Allele> alleleList, int[] baseToIndex, int[] baseCounts, int[] qualSums,
			StringBuffer baseCountsBuffer) {

		Arrays.fill(baseCounts, 0);
		int diversityCount = 0;
		ReadBackedPileup pileup = ac.getBasePileup();
		if (minBaseQ > 0) 
			pileup = pileup.getBaseFilteredPileup(minBaseQ);
		for (PileupElement e : pileup) {
			byte b = e.getBase();
			int idx = baseToIndex[b - Byte.MIN_VALUE];
			if (idx < 0)
				continue;
			if (baseCounts[idx]++ == 0)
				diversityCount++;
			qualSums[idx] += e.getQual();
		}
		List<Allele> genotypeAlleles;
		switch (genotypingMode) {
		case ONE_PER_ALLELE :
			genotypeAlleles = onePerAlleleGenotypeAlleleList(alleleList,
				baseToIndex, baseCounts, qualSums,diversityCount, baseCountsBuffer);
			break;
		case MAJOR_ALLELE :
			genotypeAlleles = majorAlleleGenotypeAlleleList(alleleList,
			    baseToIndex, baseCounts, qualSums, diversityCount, baseCountsBuffer);
			break;
		case REF_VS_ALT :
			genotypeAlleles = refVsAltAllelesGenotypeAlleleList(alleleList,
				baseToIndex, baseCounts, qualSums, diversityCount, baseCountsBuffer);
			break;
		case REF_VS_ALT_01 :
			genotypeAlleles = refVsAltAlleles01GenotypeAlleleList(alleleList,
					baseToIndex, baseCounts, qualSums, diversityCount, baseCountsBuffer);
			break;
		default:
			throw new RuntimeException("unexpected genotyping mode " + genotypingMode);
		}
		
		Map<String, Object> attributes = Collections.singletonMap(
				Constants.READ_COUNTS_KEY, (Object) baseCountsBuffer.toString());
		Genotype result = new Genotype(sampleName, genotypeAlleles,-1,null,attributes,false);
		return result;
	}

	private List<Allele> onePerAlleleGenotypeAlleleList(List<Allele> alleleList,
			int[] baseToIndex, int[] baseCounts, int[] qualSums, int diversityCount,
			StringBuffer baseCountsBuffer) {
		baseCountsBuffer.setLength(0);
		if (diversityCount == alleleList.size()) {
			for (Allele a : alleleList) {
				int baseCount = baseCounts[baseToIndex[a.getBases()[0]
						- Byte.MIN_VALUE]];
				baseCountsBuffer.append(baseCount).append(',');
			}
			baseCountsBuffer.setLength(baseCountsBuffer.length() - 1);
			return alleleList;
		} else if (diversityCount == 0) {
			return Collections.singletonList(Allele.NO_CALL);
		} else {
			Allele[] array = new Allele[diversityCount];
			int nextIndex = 0;
			for (Allele a : alleleList) {
				int baseCount = baseCounts[baseToIndex[a.getBases()[0]
						- Byte.MIN_VALUE]];
				baseCountsBuffer.append(baseCount).append(',');
				if (baseCount != 0) {
					array[nextIndex++] = a;
				}
			}
			// remove trailing comma.
			baseCountsBuffer.setLength(baseCountsBuffer.length() - 1);
			return Arrays.asList(array);
		}
	}
	
	private List<Allele> majorAlleleGenotypeAlleleList(List<Allele> alleleList,
			int[] baseToIndex, int[] baseCounts, int[] qualSums, int diversityCount,
			StringBuffer baseCountsBuffer) {
		baseCountsBuffer.setLength(0);
		int majorCount = -1;
		boolean tie = false;
		int majorQualSum = -1;
		Allele majorAllele = null;
		for (Allele a : alleleList) {
			int idx = baseToIndex[a.getBases()[0] - Byte.MIN_VALUE];
			int baseCount = baseCounts[idx];
			int qualSum = qualSums[idx];
			if (baseCount > majorCount) {
				majorAllele = a;
				majorCount = baseCount;
				majorQualSum = qualSum;
				tie = false;
			}
			else if (baseCount == majorCount) {
				switch (untieMethod) {
				case QUALITY:
				case QUALITY_AND_ORDER: 
					if (qualSum > majorQualSum) {
						majorAllele = a;
						majorCount = baseCount;
						majorQualSum = qualSum;
						tie = false;
					}
					else if (qualSum == majorQualSum) {
						tie = true;
					}
					break;
				case NO_CALL:
					tie = true;
					break;
				case ORDER:		
				}
			}
			baseCountsBuffer.append(baseCount).append(',');
		}
		baseCountsBuffer.setLength(baseCountsBuffer.length() - 1);
		if (majorAllele == null)
			return Collections.singletonList(Allele.NO_CALL);
		else if (tie)
			return Collections.singletonList(untieMethod == UntieMethod.QUALITY_AND_ORDER ? majorAllele : Allele.NO_CALL);
		else
			return Collections.singletonList(majorAllele);
	}	

	private List<Allele> refVsAltAllelesGenotypeAlleleList(List<Allele> alleleList,
			int[] baseToIndex, int[] baseCounts, int[] qualSums, int diversityCount,
			StringBuffer baseCountsBuffer) {
		baseCountsBuffer.setLength(0);
		int nonRefCount = 0;
		int refQualSum = 0;
		int nonRefQualSum = 0;
		int refCount = 0;
		int majorNonRefCount = -1;
		int majorNonRefQualSum = -1;
		Allele refAllele = null;
		Allele majorNonRefAllele = null;
		for (Allele a : alleleList) {
			int idx = baseToIndex[a.getBases()[0] - Byte.MIN_VALUE];
			int baseCount = baseCounts[idx];
			int qualSum = qualSums[idx];
			if (a.isReference()) {
				refCount += baseCount;
				refQualSum += qualSum;
				refAllele = a;
			}
			else {
				nonRefCount += baseCount;
				nonRefQualSum += qualSum;
				if (baseCount > majorNonRefCount) {
					majorNonRefAllele = a;
					majorNonRefQualSum = qualSum;
					majorNonRefCount = baseCount;
				}
				else if (baseCount == majorNonRefCount) {
					if (qualSum > majorNonRefQualSum && 
						(untieMethod == UntieMethod.QUALITY || untieMethod == UntieMethod.QUALITY_AND_ORDER)) {
						majorNonRefAllele =a;
						majorNonRefQualSum = qualSum;
						majorNonRefCount = baseCount;
					}
				}
			}
			baseCountsBuffer.append(baseCount).append(',');
		}
		baseCountsBuffer.setLength(baseCountsBuffer.length() - 1);
		if (refAllele.isSymbolic() || refAllele.isNull() || refAllele.isNoCall())
			return Collections.singletonList(majorNonRefAllele);
		if (refCount > nonRefCount)
			return Collections.singletonList(refAllele);
		else if (nonRefCount > refCount)
			return Collections.singletonList(majorNonRefAllele);
		else if (refCount == 0) {
		    return Collections.singletonList(Allele.NO_CALL);
		} 
		else {
			switch (untieMethod) {
			case NO_CALL:
				return Collections.singletonList(Allele.NO_CALL);
			case QUALITY:
			case QUALITY_AND_ORDER:
				if (refQualSum > nonRefQualSum) {
					return Collections.singletonList(refAllele);
				}
				else if (refQualSum < nonRefQualSum) {
					return Collections.singletonList(majorNonRefAllele);
				}
				else 
					return Collections.singletonList(untieMethod == UntieMethod.QUALITY ? Allele.NO_CALL : refAllele);
			case ORDER:
				return Collections.singletonList(refAllele);
			default:
				throw new RuntimeException("unexpected untie method: " + untieMethod);
			}
		}
	}	
	
	private List<Allele> refVsAltAlleles01GenotypeAlleleList(List<Allele> alleleList,
			int[] baseToIndex, int[] baseCounts, int[] qualSums, int diversityCount,
			StringBuffer baseCountsBuffer) {
		baseCountsBuffer.setLength(0);
		int nonRefCount = 0;
		int refCount = 0;
		Allele refAllele = null;
		Allele firstNonRefAllele = null;
		for (Allele a : alleleList) {
			int baseCount = baseCounts[baseToIndex[a.getBases()[0]
						- Byte.MIN_VALUE]];
			if (a.isReference()) {
				refCount += baseCount;
				refAllele = a;
			}
			else {
				nonRefCount += baseCount;
				if (firstNonRefAllele == null)
					firstNonRefAllele = a;
			}
			baseCountsBuffer.append(baseCount).append(',');
		}
		baseCountsBuffer.setLength(baseCountsBuffer.length() - 1);
		if (refAllele.isSymbolic() || refAllele.isNull() || refAllele.isNoCall())
			return Collections.singletonList(firstNonRefAllele);
		if (refCount > nonRefCount)
			return Collections.singletonList(refAllele);
		else if (nonRefCount > refCount)
			return Collections.singletonList(firstNonRefAllele);
		else 
			return Collections.singletonList(Allele.NO_CALL);
	}	

	
	
	@Override
	public boolean isReduceByInterval() {
		return false;
	}
	
	@Override
	public ReadCountStatistics reduceInit() {
		return new ReadCountStatistics();
	}

	@Override
	public ReadCountStatistics reduce(VariantContext value,
			ReadCountStatistics sum) {

		// can't call the locus because of no coverage
		if (value == null)
			return sum;

		out.add(value);
		return sum;
	}

	@Override
	public ReadCountStatistics treeReduce(ReadCountStatistics lhs,
			ReadCountStatistics rhs) {
	    // Nothing to do here for now.
		return lhs;
	}

    /**
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     *  as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).
     *  Records that are filtered in the comp track will be ignored.
     *  Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Input(fullName="comp", shortName = "comp", doc="comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }
	

}
