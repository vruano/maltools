package net.malariagen.gatk.annotators;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicReference;

import static net.malariagen.gatk.annotators.Constants.*;
import net.malariagen.gatk.coverage.LocusCategory;
import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.math.IntegerDistributionSet;
import net.malariagen.gatk.math.IntegerSampleDistributionSet;
import net.malariagen.gatk.math.IntegerSequenceDistributionSet;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class CoverageAnnotation implements InfoFieldAnnotation {

	public static AtomicReference<IntegerDistributionSet> coverageDistributionSet = new AtomicReference<IntegerDistributionSet>(
			null);

	public static Collection<String> ANNOTATION_NAMES = Collections
			.unmodifiableCollection(Arrays.asList(new String[] {
					CoverageAnnotation.class.getSimpleName(),
					CoverageAnnotation.class.getSimpleName().replace(
							"Annotation", "") }));

	private static List<String> KEY_NAMES = Collections.unmodifiableList(Arrays
			.asList(new String[] { COVERAGE_MEDIAN_FOLD_KEY,
					COVERAGE_CUMULATIVE_PROBABILITY_KEY }));

	private static List<VCFInfoHeaderLine> HEADER_LINES = Collections
			.unmodifiableList(Arrays
					.asList(new VCFInfoHeaderLine[] {
							new VCFInfoHeaderLine(COVERAGE_MEDIAN_FOLD_KEY, 1,
									VCFHeaderLineType.Float,
									"Coverage median fold factor of this locus across all samples"),
							new VCFInfoHeaderLine(
									COVERAGE_CUMULATIVE_PROBABILITY_KEY, 1,
									VCFHeaderLineType.Float,
									"Coverage cumulative probability across all samples") }));


	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
		IntegerDistributionSet cds = coverageDistributionSet.get();
		if (cds == null)
			return Collections.emptyMap();
		int depth = depth(stratifiedContexts);
		boolean isCoding = CodingAnnotation.isCoding(tracker);
		IntegerDistribution ds = resolveCoverageDistribution(ref, cds,
				isCoding);
		double cumulativeProbability = ds.cumulativeProbability(depth);
		double median = ds.median();
		double medianFoldFactor = ds.median() == 0 ? Double.NaN : ((double) depth) / median;
		Map<String, Object> map = new LinkedHashMap<String, Object>(2);
		map.put(COVERAGE_MEDIAN_FOLD_KEY, formatFractionAnnotation(medianFoldFactor));
		map.put(COVERAGE_CUMULATIVE_PROBABILITY_KEY, formatFractionAnnotation(cumulativeProbability));
		return map;
	}

	private IntegerDistribution resolveCoverageDistribution(
			ReferenceContext ref, IntegerDistributionSet cds, boolean isCoding) {
		String sequence = ref.getLocus().getLocation().getContig();
		IntegerSampleDistributionSet csds = cds.getAllSamplesDistributionSet();
		IntegerSequenceDistributionSet cseqds = csds
				.getSequenceDistributionSet(sequence);
		LocusCategory cat = isCoding ? LocusCategory.CODING
				: LocusCategory.NON_CODING;
		IntegerDistribution ds = cseqds.getCategoryDistribution(cat);
		return ds;
	}

	static int depth(Map<String, AlignmentContext> stratifiedContexts) {

		int depth = 0;
		for (Map.Entry<String, AlignmentContext> sample : stratifiedContexts
				.entrySet())
			depth += sample.getValue().size();
		return depth;
	}

	@Override
	public List<String> getKeyNames() {
		return KEY_NAMES;
	}

	@Override
	public List<VCFInfoHeaderLine> getDescriptions() {
		return HEADER_LINES;
	}

}
