package net.malariagen.gatk.coverage;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import net.malariagen.gatk.coverage.CoverageBiasWalker.GroupBy;
import net.malariagen.gatk.coverage.ReferenceComplexityWalkerWrapper.LocusComplexity;
import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.utils.ReadGroup;
import net.malariagen.gatk.utils.ReadGroupDB;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

@By(DataSource.REFERENCE)
@ReadFilters(NotPrimaryAlignmentFilter.class)
@PartitionBy(PartitionType.CONTIG)
public class CoverageBiasCorrectionWalker
		extends
		LocusWalker<net.malariagen.gatk.coverage.CoverageBiasCorrectionWalker.Locus, net.malariagen.gatk.coverage.CoverageBiasCorrectionWalker.Region> {

	private ReferenceComplexityWalkerWrapper complexity;
	private CoverageBiasCovariateCounts coverageBias;

	@Override
	public void onTraversalDone(Region result) {
		super.onTraversalDone(result);
	}

	@Override
	public boolean isReduceByInterval() {
		super.isReduceByInterval();
		return true;
	}

	@Argument(shortName = "groupBy", doc = "inidacate if the correction is per sample or read group", required = false)
	public GroupBy groupBy = GroupBy.RG;

	@Argument(shortName = "fl", fullName = "fragmentLengths", doc = "assumed fragment length or location of the directory containing fragment length distributions", required = true)
	public String fragmentLengthSpec;

	@Argument(shortName = "cb", fullName = "coverageBias", doc = "location of the directory containing coverage bias result", required = true)
	public File coverageBiasFile;

	@Argument(shortName = "W", fullName = "windowSize", required = false, doc = "size of the depth median window to correct")
	public int windowSize = 500;

	@Argument(shortName = "mmq", doc = "minimum mapping quality for a read to be considered to count depth", required = false)
	public int minMappingQuality = 0;

	@Output(shortName = "o", fullName ="ouput", doc= "VCF output file ", required = true)
	private VCFWriter writer;
	
	@Argument(shortName = "skip", doc= "number of bp between every output locus generated", required = false)
	private int skip = 1;
	
	@Argument(shortName = "burnIn", doc = "number of bp to skip at the beginning of each chromosome", required = false)
	private int burnIn = 0;
	
	
	
	private File fragmentLengthsFile;
	private Integer defaultFragmentLength;
	private FragmentLengthSummary fragmentLengths;
	private CoverageBiasCovariateCounts coverageBiasCounter;

	private Set<String> groupNames;
	private Map<String, Integer> groupIndex;
	private HashMap<String, Integer> fragmentSizeByName;
	private ReadGroupDB readGroupDb;

	@Override
	public void initialize() {
		super.initialize();
		if (groupBy == GroupBy.NONE || groupBy == GroupBy.WS) {
			throw new IllegalArgumentException(
					"only SMRG, RG or SM allowed for the groupBy option");
		}
		if (skip <= 0)
			throw new UserException("skip cannot 0 or negative");
		
		if (burnIn < 0)
			throw new UserException("burnIn cannot be negative");
		
		readGroupDb = new ReadGroupDB(getToolkit());
		groupNames = groupBy.buildGroupNameSet(readGroupDb);
		groupIndex = new HashMap<String,Integer>(groupNames.size());
		int nextIdx = 0;
		for (String gn : groupNames)
			groupIndex.put(gn,nextIdx++);
		initializeCoverageBiasCounter();
		initializeFragmentLengths();

		try {
			coverageBias = CoverageBiasCovariateCounts
					.loadFrom(coverageBiasFile);
		} catch (IOException e) {
			throw new StingException(
					"problems loading coverage bias counts file", e);
		}

		complexity = new ReferenceComplexityWalkerWrapper(getToolkit(),
				groupBy, fragmentLengthsFile, defaultFragmentLength,
				coverageBias.getGroupNames());
		
		Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		headerLines.add(new VCFFormatHeaderLine("FFS",1,VCFHeaderLineType.Integer,""));
		headerLines.add(new VCFFormatHeaderLine("FFSr",1,VCFHeaderLineType.Float,""));
		headerLines.add(new VCFFormatHeaderLine("FFSc",1,VCFHeaderLineType.Float,""));
		headerLines.add(new VCFFormatHeaderLine("FFSx",1,VCFHeaderLineType.Float,""));
		headerLines.add(new VCFFormatHeaderLine("FGC",1,VCFHeaderLineType.Float,""));
		headerLines.add(new VCFFormatHeaderLine("WS",1,VCFHeaderLineType.Integer,""));
				VCFHeader header = new VCFHeader(headerLines,groupNames);

		writer.writeHeader(header);
	}

	private void initializeCoverageBiasCounter() {
		if (!coverageBiasFile.exists())
			throw new UserException(
					"you must provide a existing coverage bias analysis output file location");
		try {
			coverageBiasCounter = CoverageBiasCovariateCounts
					.loadFrom(coverageBiasFile);
		} catch (IOException e) {
			throw new StingException(e.getMessage(), e);
		}
	}

	private void initializeFragmentLengths() {

		File candidateFile = new File(this.fragmentLengthSpec);
		if (candidateFile.exists()) {
			fragmentLengthsFile = candidateFile;
			try {
				fragmentLengths = FragmentLengthSummary
						.loadFrom(fragmentLengthsFile);
			} catch (IOException e) {
				throw new StingException(e.getMessage(), e);
			}
		} else if (fragmentLengthSpec.matches("(+|-)\\s*\\d+")) {
			try {
				defaultFragmentLength = Integer.parseInt(fragmentLengthSpec);
				if (defaultFragmentLength <= 0)
					throw new UserException(
							"the fragment length provided must be greater than 0 ("
									+ fragmentLengthSpec + ")");
			} catch (NumberFormatException e) {
				throw new UserException(
						"the fragment length provided does not seem to be a valid integer number ("
								+ fragmentLengthSpec + ")");
			}
		} else {
			throw new UserException(
					"you must indicate a fragment length fixed size or a distribution file/directory location (-fl option)");
		}

		fragmentSizeByName = new HashMap<String, Integer>();
		for (String gn : groupNames) {
			if (fragmentLengths != null) {
				ReadGroup rg = readGroupDb.getReadGroupByID(gn);
				IntegerDistribution id = null;
				if (rg != null)
					id = fragmentLengths.getReadGroupFragmentLengths(gn);
				else {
					Sample s = readGroupDb.getSampleDB().getSample(gn);
					if (s != null)
						id = fragmentLengths
								.getSampleFragmentLengths(s.getID());
				}
				int size = (id == null) ? (defaultFragmentLength == null ? -1
						: defaultFragmentLength) : (int) Math
						.round(id.median());
				fragmentSizeByName.put(gn, size);
			}
		}
	}

	@Override
	public Locus map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		Locus result = new Locus();
		result.refContext = ref;
		complexity.map(tracker, ref, context);
		Map<String, AlignmentContext> stratifiedContexts = CoverageBiasWalker
				.stratifyByGroupName(context, groupNames, groupBy);
		result.addCounts(stratifiedContexts);
		return result;
	}

	@Override
	public Region reduceInit() {
		Region result = new Region();
		result.complexity = complexity.reduceInit();
		return result;
	}

	@Override
	public Region reduce(Locus value, Region sum) {
		sum.add(value);
	
		sum.complexity = complexity.reduce(value.refContext, sum.complexity);
		List<LocusComplexity> completed = complexity.removeCompleted();
		GenomeLoc lastLoc = null;
		for (LocusComplexity lc : completed) {
			Locus l = sum.get(lc.locus);
			l.setComplexity(lc);
			lastLoc = lc.locus;
		}
		if (lastLoc != null)
			sum.emitTo(lastLoc);
		return sum;
	}

	public class Region {
		public MultiWindowSequenceComplexity complexity;

		public SortedMap<GenomeLoc, Locus> locusBuffer = new TreeMap<GenomeLoc, Locus>();


		public void add(Locus value) {
			if (locusBuffer.size() > 0 && locusBuffer.firstKey().getContig() != value.refContext.getLocus().getContig())
				locusBuffer.clear();
			locusBuffer.put(value.getLocation(), value);
		}

		public void emitTo(GenomeLoc toLoc) {
			if (locusBuffer.isEmpty())
				return;
			GenomeLoc lastLoc = null;
			for (Locus l : locusBuffer.values()) {
				GenomeLoc loc = l.getLocation();
				lastLoc = loc;
				if (toLoc.getStart() - loc.getStart() + 1 < windowSize)
					break;
				Iterator<Locus> it = locusBuffer.tailMap(loc).values()
						.iterator();
				if (burnIn > 0 && loc.getStart() <= burnIn)
					continue;
				if (skip > 1 && ((loc.getStart() - burnIn) % skip) != 1)
					continue;
				for (int i = 0; i < windowSize; i++) {
					Locus m = it.next();
					for (int j = 0; j < l.byGroupIndex.length; j++) {
						if (Double.isNaN(m.byGroupIndex[j].windowForwardLambda) || Double.isNaN(m.byGroupIndex[j].forwardLambda))
							continue;
					//	if (m.byGroupIndex[j].depth < 5)
					//		continue;
						l.byGroupIndex[j].windowSize++;
						l.byGroupIndex[j].windowForwardStartsCorrected += Math.max(m.byGroupIndex[j].forwardStarts,0.1) / m.byGroupIndex[j].forwardLambda;
						l.byGroupIndex[j].windowForwardGcBias += m.byGroupIndex[j].forwardGcBias;
						l.byGroupIndex[j].windowReverseStarts += m.byGroupIndex[j].reverseStarts;
						l.byGroupIndex[j].windowForwardStarts += m.byGroupIndex[j].forwardStarts;
						l.byGroupIndex[j].windowForwardLambda += m.byGroupIndex[j].forwardLambda;
						l.byGroupIndex[j].windowReverseCorrected += m.byGroupIndex[j].reverseStartCorrected;
					}
				}
				VariantContext vc = l.toVariantContext();
				writer.add(vc);
			}
			locusBuffer.headMap(lastLoc).clear();

		}

		public Locus get(GenomeLoc loc) {
			Locus l = locusBuffer.get(loc);
			if (l == null)
				throw new IllegalStateException(
						"locus complexity reported before locus has been visited!!!");
			return l;
		}

	}

	public class LocusGenotype {

		public int windowSize;
		public int forwardStarts;
		public int reverseStarts;
		public int depth;
		public double forwardLambda;
		public double reverseStartCorrected;
		public int windowForwardStarts;
		public int windowReverseStarts;
		public int windowForwardStartsCorrected;
		public double windowForwardLambda;
		public double windowReverseCorrected;
		public double forwardGcBias;
		public double reverseGcBias;
		public double windowForwardGcBias;

	}

	public class Locus {

		public ReferenceContext refContext;
		LocusGenotype[] byGroupIndex = new LocusGenotype[groupNames.size()];

		public GenomeLoc getLocation() {
			if (refContext == null)
				return null;
			return refContext.getLocus();
		}

		public VariantContext toVariantContext() {
			VariantContextBuilder vcb = new VariantContextBuilder();
			vcb.loc(refContext.getLocus());
			vcb.alleles(Collections.singletonList(Allele.create(refContext.getBase(),true)));
			GenotypesContext gc = GenotypesContext.create();
			List<Allele> noCall = Collections.emptyList();
			Set<String> noFilters = Collections.emptySet();
			for (String gn : groupNames) {
				int idx = groupIndex.get(gn);
				Map<String,Object> attributes = new LinkedHashMap<String,Object>(10);
				LocusGenotype lg = byGroupIndex[idx];
				if (lg == null) {
					gc.add(new Genotype(gn,noCall));
					continue;
				}
				attributes.put("WS",lg.windowSize);
				attributes.put("FGC",lg.windowForwardGcBias / (double) lg.windowSize);
				attributes.put("RGC",lg.reverseGcBias);
				attributes.put("FFS",2.0 * ((double) lg.windowForwardStarts) / (double) lg.windowSize);
				attributes.put("FFSr",lg.windowForwardLambda / (double) lg.windowSize);
				attributes.put("FFSc",2.0 * ((double)lg.windowForwardStarts) / (double) lg.windowForwardLambda);
				attributes.put("FFSx",2.0 * ((double)lg.windowForwardStartsCorrected) / (double) lg.windowSize);
				//attributes.put("RFS",lg.windowReverseStarts / (double) lg.windowSize);
				//attributes.put("RFSr",lg.windowReverseCorrected / (double) lg.windowSize);
				//attributes.put("RFSc",((double)lg.windowReverseStarts) / lg.windowReverseCorrected);
				//attributes.put("BFS", ((double) lg.windowForwardStarts + lg.windowReverseStarts ) / (double) lg.windowSize);
				//attributes.put("BFSc",(((double) lg.windowForwardStarts)  / lg.windowForwardCorrected 
				//		  + ((double) lg.windowReverseStarts) / lg.windowReverseCorrected));
				gc.add(new Genotype(gn,noCall,1.0,noFilters,attributes,false));
			}
			vcb.genotypes(gc);
			return vcb.make();
		}

		public void addCounts(Map<String, AlignmentContext> stratifiedContexts) {
			for (Map.Entry<String, AlignmentContext> e : stratifiedContexts
					.entrySet()) {
				AlignmentContext ac = e.getValue();
				String name = e.getKey();
				int index = groupIndex.get(name);
				LocusGenotype lg = byGroupIndex[index];
				if (lg == null)
					byGroupIndex[index] = lg = new LocusGenotype();
				int forwardStarts = 0;
				int reverseStarts = 0;
				for (PileupElement pe : ac.getBasePileup()) {
					if (pe.isDeletion() || pe.isInsertionAtBeginningOfRead())
						continue;
					lg.depth++;
					if (pe.getMappingQual() < minMappingQuality)
						continue;
					SAMRecord read = pe.getRead();
					if (!read.getFirstOfPairFlag())
						continue;
					if (pe.getOffset() == 0) {
						if (!read.getReadNegativeStrandFlag())
							forwardStarts++;
					} else if (pe.getOffset() == pe.getRead().getReadLength() - 1) {
						if (read.getReadNegativeStrandFlag())
							reverseStarts++;
					}
				}
				lg.forwardStarts += forwardStarts;
				lg.reverseStarts += reverseStarts;
			}
		}

		public void setComplexity(LocusComplexity lc) {
			for (String gn : groupNames) {
				int i = groupIndex.get(gn);
				LocusGenotype lg = byGroupIndex[i];
				lg.forwardGcBias = lc.getGcBias(gn, true);
				lg.reverseGcBias = lc.getGcBias(gn, false);
				lg.forwardLambda = coverageBiasCounter.getStartRate(gn,
								fragmentSizeByName.get(gn) -1,
								lc.getGcBias(gn, true),true);
				lg.reverseStartCorrected = coverageBiasCounter.getStartRate(gn,
								fragmentSizeByName.get(gn) -1,
								lc.getGcBias(gn, false),false);
			}

		}
	}

}
