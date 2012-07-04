package net.malariagen.gatk.coverage;

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
import net.malariagen.gatk.coverage.SequenceComplexity.LocusComplexity;
import net.malariagen.vcf.VCFReadGroupHeaderLine;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.utils.GenomeLoc;
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
public class ReferenceComplexityWalker extends
		LocusWalker<ReferenceContext, MultiWindowSequenceComplexity> {

	@Argument(shortName = "W", fullName = "windowSize", required = false, doc = "window-size to get stats on")
	protected Integer windowSize = null;

	@Argument(shortName = "fs", fullName = "fragmentLenghts", required = false, doc = "fragment-length output for the samples under analysis. A genotype column will be output with the corresponding median window size")
	protected File fragmentLengthsFile = null;

	@Argument(shortName = "rounding", fullName = "windowSizeRounding", doc="window size rounding, 1 means no rounding", required=false)
	protected int rounding = 1;
	
	@Output(shortName = "o", fullName = "output", doc = "File to which variants should be written", required = true)
	protected VCFWriter writer = null;

	
	@Override
	public ReferenceContext map(RefMetaDataTracker tracker,
			ReferenceContext ref, AlignmentContext context) {
		return ref;
	}

	@Output(shortName = "groupBy", required = false)
	protected GroupBy groupBy = null;

	private FragmentLengthSummary fragmentLengthSummary;

	private Map<String, Integer> groupWindowSize;

	@Override
	public void initialize() {
		super.initialize();
		if (writer == null)
			throw new IllegalStateException("the write is yet null");
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

		Set<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		VCFHeader header = null;
		if (groupBy == null)
			groupBy = fragmentLengthsFile == null ? GroupBy.NONE : GroupBy.SMRG;
		if (groupBy != GroupBy.NONE && fragmentLengthsFile == null)
			throw new UserException(
					"if you request a grouping (groupBy != NONE) you need also to provide a fragment-lenghts location alwell (-fs LOC)");
		if (groupBy == GroupBy.NONE && windowSize == null)
			throw new UserException(
					"if you are not grouping output by sample you must specify a common window size (-W NUM)");

		groupWindowSize = new HashMap<String, Integer>(100);
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
			headerLines
					.add(new VCFFormatHeaderLine("ED", 1,
							VCFHeaderLineType.Float,
							"Stop positio of the interval considerde in complex INFO fields"));

			Set<String> groupNames = new LinkedHashSet<String>(100);
			if (groupBy.implies(GroupBy.RG) || groupBy == groupBy.WS) {
				groupNames.addAll(fragmentLengthSummary.getReadGroups());
				for (String name : groupNames) {
					int ws = (int) Math.round(fragmentLengthSummary
							.getReadGroupFragmentLengths(name).median() / rounding ) * rounding;
					groupWindowSize.put(name, ws);
					headerLines.add(new VCFReadGroupHeaderLine(name,
							Collections.singletonMap("WindowSize", "" + ws)));
				}
			}
			if (groupBy.implies(GroupBy.SM) || groupBy == GroupBy.WS) {
				groupNames.addAll(fragmentLengthSummary.getSamples());
				for (String name : fragmentLengthSummary.getSamples())
					groupWindowSize.put(name, (int) Math
							.round(fragmentLengthSummary
									.getSampleFragmentLengths(name).median() / rounding) * rounding);
			}
			if (groupBy == GroupBy.WS) {
				Set<String> realGroupNames = groupNames;
				groupNames = new HashSet<String>(realGroupNames.size());
				for (String name : realGroupNames) {
					int ws = groupWindowSize.remove(name);
					groupNames.add("" + ws);
					groupWindowSize.put("" + ws,ws);
				}
			}
			header = new VCFHeader(headerLines, groupNames);
		}

		if (windowSize != null) {
			headerLines.add(new VCFInfoHeaderLine("GCBias", 1,
					VCFHeaderLineType.Float,
					"GC bias expressed in a percentage from 0 to 100"));
			headerLines.add(new VCFInfoHeaderLine("NucEnt", 1,
					VCFHeaderLineType.Float,
					"Nucleotide entropy using Euler's e as base"));
			headerLines.add(new VCFInfoHeaderLine("TriEnt", 1,
					VCFHeaderLineType.Float,
					"Trinucleotide entropy using Euler's e as base"));
			headerLines
					.add(new VCFInfoHeaderLine("END", 1,
							VCFHeaderLineType.Float,
							"Stop positio of the interval considerde in complex INFO fields"));
			header = new VCFHeader(headerLines);
		}
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
		return new MultiWindowSequenceComplexity(windowSizeInts);
	}

	@Override
	public MultiWindowSequenceComplexity reduce(ReferenceContext value,
			MultiWindowSequenceComplexity sum) {
		List<Map<Integer, SequenceComplexity.LocusComplexity>> lcm = sum
				.count(value);
		for (Map<Integer, SequenceComplexity.LocusComplexity> l : lcm)
			if (l.size() != 0)
				emit(l);
		return sum;
	}

	private void emit(Map<Integer, LocusComplexity> lcm) {
		VariantContextBuilder vcb = new VariantContextBuilder();
		LocusComplexity example = lcm.values().iterator().next();
		Collection<Allele> alleles = Collections.singleton(Allele.create(
				example.getRefNuc().byteValue(), true));
		vcb.alleles(alleles);
		GenomeLoc loc = example.getLocus();
		if (windowSize != null && lcm.containsKey(windowSize)) {
			LocusComplexity lc = lcm.get(windowSize);
			Map<String, Object> attributes = new LinkedHashMap<String, Object>(
					4);
			attributes.put("GCBias", lc.getGcBias());
			attributes.put("NucEnt", lc.getNucEnt());
			attributes.put("TriEnt", lc.getTriEnt());
			attributes.put("END", lc.getLocus().getStart() + lc.size() - 1);
			vcb.attributes(attributes);
			if (lc.getLocus().compareTo(loc) < 0)
				loc = lc.getLocus();
		}
		List<Allele> noCall = Collections.singletonList(Allele.NO_CALL);
		Set<String> noFilters = Collections.emptySet();
		GenotypesContext gc = GenotypesContext.create();
		for (Map.Entry<String, Integer> gws : groupWindowSize.entrySet()) {
			LocusComplexity lc = lcm.get(gws.getValue());
			Map<String, Object> attr = new LinkedHashMap<String, Object>(4);
			if (lc != null) {
				attr.put("GC", lc.getGcBias());
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

}
