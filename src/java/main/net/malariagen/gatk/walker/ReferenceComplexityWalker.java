package net.malariagen.gatk.walker;

import java.io.File;
import java.io.IOException;
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
import net.malariagen.gatk.utils.ReadGroup;
import net.malariagen.gatk.utils.ReadGroupDB;
import net.malariagen.gatk.walker.SequenceComplexity.LocusComplexity;
import net.malariagen.vcf.VCFReadGroupHeaderLine;
import net.sf.samtools.SAMReadGroupRecord;

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
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
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

	static class Locus {
		ReferenceContext ref;
		int refMQ;
	}

	@Argument(shortName = "W", fullName = "windowSize", required = false, doc = "window-size to get stats on")
	protected List<Integer> windowSizeList = Collections.emptyList();

	@Argument(shortName = "useInfoFields", fullName = "useInfoFields", required = false, doc = "whether complexity for a fixed unique window size shoul be outputted in the info fields")
	protected boolean useInfoFields = false;

	protected Integer windowSize;

	@Argument(shortName = "fs", fullName = "fragmentLenghts", required = false, doc = "fragment-length output for the samples under analysis. A genotype column will be output with the corresponding median window size")
	protected File fragmentLengthsFile = null;

	@Argument(shortName = "rounding", fullName = "windowSizeRounding", doc = "window size rounding, 1 means no rounding", required = false)
	protected int rounding = 1;

	@Output(shortName = "o", fullName = "output", doc = "File to which variants should be written", required = true)
	protected VCFWriter writer = null;

	@Argument(shortName = "exaustiveRefReads", fullName = "exaustiveReferenceReads", doc = "name of read-group or sample containing exaustive reference reads", required = false)
	protected String exaustiveRefReadGroupOrSample = null;

	@Override
	public Locus map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		Locus result = new Locus();
		result.ref = ref;
		if (exaustiveRef != null) {
			for (PileupElement pe : context.getBasePileup()) {
				GATKSAMRecord r = pe.getRead();
				SAMReadGroupRecord rg = r.getReadGroup();
				if (r.getAlignmentStart() != ref.getLocus().getStart())
					continue;
				if (!rg.getId().equals(this.exaustiveRefReadGroupOrSample)
						&& !rg.getSample().equals(
								this.exaustiveRefReadGroupOrSample))
					continue;
				String contig = ref.getLocus().getContig();
				String readName = r.getReadName();
				if (!readName.startsWith(contig))
					continue;
				int idx = contig.length();
				StringBuffer sb = new StringBuffer(10);
				boolean numberFound = false;
				while (idx < readName.length()) {
					char c = readName.charAt(idx++);
					if (Character.isDigit(c)) {
						numberFound = true;
						sb.append(c);
					} else if (numberFound)
						break;
				}
				if (!numberFound)
					continue;
				int num = Integer.parseInt(sb.toString());
				if (num == ref.getLocus().getStart()) {
					result.refMQ = r.getMappingQuality();
					break;
				}
			}
		}
		return result;
	}

	@Output(shortName = "groupBy", required = false)
	protected GroupBy groupBy = null;

	private FragmentLengthSummary fragmentLengthSummary;

	private ReadGroupDB readGroupDb;

	Map<String, Integer> groupWindowSize;
	Map<String, Integer> realGroupWindowSize;

	private Integer exaustiveRef;

	@Override
	public void initialize() {
		super.initialize();
		if (writer == null)
			throw new IllegalStateException("the write is yet null");
		if (!windowSizeList.isEmpty())
			windowSize = windowSizeList.get(0);
		if (windowSize != null && windowSize <= 0)
			throw new IllegalArgumentException(
					"a strictly positive window size must be indicated, the one provided was "
							+ windowSize + ")");
		if (fragmentLengthsFile != null) {
			if (!fragmentLengthsFile.exists())
				throw new UserException(
						"the fragment-lengths file/directory provided does not exists: '"
								+ fragmentLengthsFile + "'");
			if (!fragmentLengthsFile.canRead())
				throw new UserException(
						"the fragment-lengths file/directory provided does not exists: '"
								+ fragmentLengthsFile + "'");
			try {
				fragmentLengthSummary = FragmentLengthSummary
						.loadFrom(fragmentLengthsFile);
			} catch (IOException e) {
				throw new UserException(
						"could not load the fragment lengths from the location provided: '"
								+ fragmentLengthsFile + "'", e);
			}
		}
		if (rounding <= 0)
			throw new UserException("the window size rounding must be positive");

		readGroupDb = new ReadGroupDB(this.getToolkit());

		if (useInfoFields && windowSizeList.size() != 1)
			throw new UserException(
					"you cannot request complexity to be output in info fields unless you indicate one and only one window size with -W");

		if (fragmentLengthsFile == null && groupBy == null)
			groupBy = GroupBy.WS;

		Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		VCFHeader header = null;
		if (groupBy == null)
			groupBy = fragmentLengthsFile == null ? GroupBy.NONE : GroupBy.SMRG;
//		if (groupBy != GroupBy.NONE && fragmentLengthsFile == null)
//			throw new UserException(
//					"if you request a grouping (groupBy != NONE) you need also to provide a fragment-lenghts location alwell (-fs LOC)");
		if (groupBy == GroupBy.NONE && windowSize == null)
			throw new UserException(
					"if you are not grouping output by sample you must specify a common window size (-W NUM)");

		if (groupBy == GroupBy.RG || GroupBy.SM == groupBy)
			if (fragmentLengthsFile == null) throw new UserException("if you request grouping per read group or sample you need to provide a fragments -lenght file");
		
		groupWindowSize = new HashMap<String, Integer>(100);
		Set<String> groupNames = new LinkedHashSet<String>(100);
		if (groupBy != GroupBy.NONE) {
			headerLines.add(new VCFFormatHeaderLine("GC", 1,
					VCFHeaderLineType.Float,
					"GC bias expressed in a percentage from 0 to 100"));
			headerLines.add(new VCFFormatHeaderLine("NE", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy using Euler's e as base"));
			headerLines.add(new VCFFormatHeaderLine("TE", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy using Euler's e as base"));
			headerLines.add(new VCFFormatHeaderLine("NC", 1,
					VCFHeaderLineType.Integer,
					"Nucleotide count on the windows (excludes 'N' 'X')"));
			headerLines
					.add(new VCFFormatHeaderLine("ED", 1,
							VCFHeaderLineType.Float,
							"Stop positio of the interval considerde in complex INFO fields"));

			if (groupBy.implies(GroupBy.RG) && fragmentLengthSummary != null) {
				groupNames.addAll(fragmentLengthSummary.getReadGroups());
				for (String name : groupNames) {
					int ws = (int) Math.round(fragmentLengthSummary
							.getReadGroupFragmentLengths(name).median()
							/ rounding) * rounding;
					groupWindowSize.put(name, ws);
					headerLines.add(new VCFReadGroupHeaderLine(name,
							Collections.singletonMap("WindowSize", "" + ws)));
				}
			}
			if (groupBy.implies(GroupBy.SM) && fragmentLengthSummary != null) {
				groupNames.addAll(fragmentLengthSummary.getSamples());
				for (String name : fragmentLengthSummary.getSamples())
					groupWindowSize.put(
							name,
							(int) Math.round(fragmentLengthSummary
									.getSampleFragmentLengths(name).median()
									/ rounding) * rounding);
			}
			Set<String> realGroupNames = groupNames;
			realGroupWindowSize = groupWindowSize;
			if (groupBy == GroupBy.WS) {
				groupNames = new LinkedHashSet<String>(realGroupNames.size());
				groupWindowSize = new LinkedHashMap<String, Integer>(
						realGroupWindowSize.size());
				for (String name : realGroupNames) {
					int ws = realGroupWindowSize.get(name);
					groupNames.add("" + ws);
					groupWindowSize.put("" + ws, ws);
				}
				if (windowSizeList.size() > 0)
					for (Integer ws : windowSizeList.subList(useInfoFields ? 1
							: 0, windowSizeList.size())) {
						groupNames.add("" + ws);
						groupWindowSize.put("" + ws, ws);
					}
			}
		}
		initializeExaustiveReference();
		if (exaustiveRef != null)
			headerLines.add(new VCFInfoHeaderLine("RefMQ", 1,
					VCFHeaderLineType.Integer,
					"Mapping quality of the exact reference read"));
		if (useInfoFields) {
			headerLines.add(new VCFInfoHeaderLine("GCBias", 1,
					VCFHeaderLineType.Float,
					"GC bias expressed in a percentage from 0 to 100"));
			headerLines.add(new VCFInfoHeaderLine("NucEnt", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy using Euler's e as base"));
			headerLines.add(new VCFInfoHeaderLine("TriEnt", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy using Euler's e as base"));
			headerLines.add(new VCFInfoHeaderLine("NucCnt", 1,
					VCFHeaderLineType.Integer,
					"Number of nucleotides in the window"));
			headerLines
					.add(new VCFInfoHeaderLine("END", 1,
							VCFHeaderLineType.Float,
							"Stop positio of the interval considerde in complex INFO fields"));
		}
		header = new VCFHeader(headerLines,groupNames);
		writer.writeHeader(header);
	}

	private void initializeExaustiveReference() {
		if (exaustiveRefReadGroupOrSample != null) {
			Set<ReadGroup> exaustiveRefReadGroups = readGroupDb
					.getReadGroupsBySampleOrReadGroupID(exaustiveRefReadGroupOrSample);
			if (exaustiveRefReadGroups.size() == 0)
				throw new UserException(
						"The exaustive reference sample or read-group provided is not found amongst the input data read-groups");
			else if (exaustiveRefReadGroups.size() == 1)
				exaustiveRef = groupWindowSize.get(exaustiveRefReadGroups
						.iterator().next());
			else {
				for (ReadGroup rg : exaustiveRefReadGroups) {
					Integer candidate = groupWindowSize.get(rg.getID());
					if (candidate == null)
						;
					else if (exaustiveRef == null)
						exaustiveRef = candidate;
					else if (exaustiveRef.intValue() != candidate.intValue())
						throw new UserException(
								"the exaustive reference sample/ read-group maps to more than one windows-size");
				}
			}
			if (exaustiveRef == null)
				exaustiveRef = windowSize;
			if (exaustiveRef == null)
				throw new UserException(
						"cannot resolve the window-size for the exaustive reference sample or read-group");
		}
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
		return new MultiWindowSequenceComplexity(windowSizeInts);
	}

	@Override
	public MultiWindowSequenceComplexity reduce(Locus value,
			MultiWindowSequenceComplexity sum) {
		if (sum.lastLocus() != null & value.ref.getLocus().compareContigs(sum.lastLocus()) != 0) 
			for (Map<Integer, SequenceComplexity.LocusComplexity> l : sum.flush())
				if (l.size() != 0) emit(l);
		List<Map<Integer, SequenceComplexity.LocusComplexity>> lcm = sum.count(
				value.ref, exaustiveRef, value.refMQ);
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
			attributes.put("GCBias", lc.getGcBias());
			attributes.put("NucEnt", lc.getNucEnt());
			attributes.put("TriEnt", lc.getTriEnt());
			attributes.put("NucCnt", lc.getNucCount());
			attributes.put("END", lc.getLocus().getStart() + lc.size() - 1);
			if (lc.getLocus().compareTo(loc) < 0)
				loc = lc.getLocus();
		}
		if (exaustiveRef != null) {
			attributes.put("RefMQ", lcm.get(exaustiveRef).getRefMQ());
		}
		vcb.attributes(attributes);
		List<Allele> noCall = Collections.singletonList(Allele.NO_CALL);
		Set<String> noFilters = Collections.emptySet();
		GenotypesContext gc = GenotypesContext.create();
		for (Map.Entry<String, Integer> gws : groupWindowSize.entrySet()) {
			LocusComplexity lc = lcm.get(gws.getValue());
			Map<String, Object> attr = new LinkedHashMap<String, Object>(4);
			if (lc != null) {
				attr.put("GC", lc.getGcBias());
				attr.put("NC", lc.getNucCount());
				attr.put("NE", lc.getNucEnt());
				attr.put("TE", lc.getTriEnt());
				attr.put("ED", lc.getLocus().getStart() + lc.size() - 1);
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
				rc, null, 0);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : lcm)
			if (l.size() != 0)
				emit(l);
		return sum;
	}

}
