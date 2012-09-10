package net.malariagen.gatk.coverage;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static net.malariagen.gatk.annotators.Constants.*;
import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.math.IntegerSampleDistributionSet;
import net.malariagen.gatk.math.IntegerSequenceDistributionSet;

import org.apache.commons.math.stat.descriptive.rank.Median;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;

public class AccesibilityAnnotations {

	private static final int MAX_MAPPING_QUALITY = 256;

	public static Collection<String> ANNOTATION_NAMES = Collections
			.unmodifiableCollection(Arrays.asList(new String[] {
					AccesibilityAnnotations.class.getSimpleName(),
					AccesibilityAnnotations.class.getSimpleName().replace(
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

	private AccessibilityWalker walker;

	private int maxAllMappingQuality;
	private int[] allMappingQualityCounts;
	private int[][] mappingQualityCounts;
	private int[] maxMappingQuality;
	private int[] groupDepth;
	private int allDepth;
	
	private Median medianCalculator = new Median();

	private Map<String, Integer> groupIndex;

	private double[] groupDoubleAuxiliar;
	private double[] groupDoubleAuxiliar2;

	public AccesibilityAnnotations(AccessibilityWalker walker) {
		if (walker == null)
			throw new IllegalArgumentException();
		this.walker = walker;
		int nextIdx = 1;
		for (String gn : walker.groupNames)
			groupIndex.put(gn, nextIdx++);
		mappingQualityCounts = new int[groupIndex.size()][MAX_MAPPING_QUALITY + 1];
		maxMappingQuality = new int[groupIndex.size()];
		maxAllMappingQuality = 0;
		allMappingQualityCounts = new int[MAX_MAPPING_QUALITY + 1];
		groupDepth = new int[groupIndex.size()];
		groupDoubleAuxiliar = new double[groupIndex.size()];
		groupDoubleAuxiliar2 = new double[groupIndex.size()];
		allDepth = 0;
	}

	public void annotate(RefMetaDataTracker tracker, ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedAlignmentContext, Map<String,Object> infoAnnotations, Map<String,Map<String,Object>> formatAnnotations) {

		// clear counters.
		Arrays.fill(allMappingQualityCounts, 0, maxAllMappingQuality, 0);
		maxAllMappingQuality = 0;
		for (int i = 0; i < groupIndex.size(); i++)
			Arrays.fill(mappingQualityCounts[i], 0, maxMappingQuality[i], 0);
		Arrays.fill(maxMappingQuality, 0);
		Arrays.fill(groupDepth, 0);
		allDepth = 0;

		for (Map.Entry<String, AlignmentContext> gn : stratifiedAlignmentContext
				.entrySet()) {
			int index = groupIndex.get(gn.getKey());
			int[] mappingQualityCount = mappingQualityCounts[index];
			allDepth += groupDepth[index] = gn.getValue().getBasePileup()
					.depthOfCoverage();
			for (PileupElement pe : gn.getValue().getBasePileup()) {
				int mq = pe.getMappingQual();
				if (mappingQualityCount[mq]++ == 0
						&& mq > maxMappingQuality[index])
					maxMappingQuality[index] = mq;
				if (allMappingQualityCounts[mq]++ == 0
						&& mq > maxAllMappingQuality)
					maxAllMappingQuality = mq;
			}
		}

		infoAnnotations.put("DP", allDepth);
		infoAnnotations.put("MQ0", allMappingQualityCounts[0]);
		infoAnnotations.put("MQ", rms(allMappingQualityCounts,0,maxAllMappingQuality + 1));
		
		for (int i = 0; i < groupDepth.length; i++)
			groupDoubleAuxiliar[i] = (double) groupDepth[i];
		
		double groupMedianDepth = medianCalculator.evaluate(groupDoubleAuxiliar);
		infoAnnotations.put("MedDP",groupMedianDepth);
		
		
		
		boolean isCoding = false;
		if (walker.features != null) {
			isCoding = walker.isCoding(tracker);
			if (isCoding)
				infoAnnotations.put("CODING", null);
		}
		if (walker.normalize) {
			IntegerSampleDistributionSet allSet = walker.coverageDistributionSet
					.getAllSamplesDistributionSet();
			IntegerSequenceDistributionSet allSeqSet = walker.normalizePerChromosome ? allSet
					.getSequenceDistributionSet(ref.getLocus().getContig())
					: allSet.getAllSequencesDistributionSet();
			IntegerDistribution allDist = walker.normalizePerCodingStatus ? allSeqSet
					.getCategoryDistribution(isCoding ? LocusCategory.CODING
							: LocusCategory.NON_CODING) : allSeqSet
					.getAllCategoriesDistribution();

			double allMedian = allDist.median();
			double allCcp = allDist.cumulativeProbability(allDepth);
			infoAnnotations.put("Cmf", allDepth / allMedian);
			infoAnnotations.put("Ccp", allCcp);
			for (String gn : walker.groupNames) {
				IntegerSampleDistributionSet gSet = walker.coverageDistributionSet
						.getSampleDistributionSet(gn);
				IntegerSequenceDistributionSet gSeqSet = walker.normalizePerChromosome ? gSet
						.getSequenceDistributionSet(ref.getLocus().getContig())
						: gSet.getAllSequencesDistributionSet();
				IntegerDistribution gDist = walker.normalizePerCodingStatus ? gSeqSet
						.getCategoryDistribution(isCoding ? LocusCategory.CODING
								: LocusCategory.NON_CODING) : gSeqSet
						.getAllCategoriesDistribution();
				int index = groupIndex.get(gn);
				Map<String,Object> annotations = formatAnnotations.get(gn);
				annotations.put("MQ0", mappingQualityCounts[index][0]);
				groupDoubleAuxiliar[index] = 100.0 * mappingQualityCounts[index][0] / groupDoubleAuxiliar[index];
				annotations.put("MQ", rms(mappingQualityCounts[index],0,maxMappingQuality[index] + 1));
				annotations.put("DP", groupDepth[index]);
				double cmf = groupDepth[index] / gDist.median();
				annotations.put("CMF", cmf);
				groupDoubleAuxiliar2[index] = cmf;
				annotations.put("CCP", gDist.cumulativeProbability(groupDepth[index]));
			}	
			double cmfMedian = medianCalculator.evaluate(groupDoubleAuxiliar2);
			infoAnnotations.put("MedCMF",cmfMedian);
		}
		else {
			for (String gn : walker.groupNames) {
				int index = groupIndex.get(gn);
				Map<String,Object> annotations = formatAnnotations.get(gn);
				annotations.putAll(annotations);
				annotations.put("MQ0", mappingQualityCounts[index][0]);
				groupDoubleAuxiliar[index] = 100.0 * mappingQualityCounts[index][0] / groupDoubleAuxiliar[index];
				annotations.put("MQ", rms(mappingQualityCounts[index],0,maxMappingQuality[index] + 1));
				annotations.put("DP", groupDepth[index]);
			}	
		}
		double mq0Median = medianCalculator.evaluate(groupDoubleAuxiliar);
		infoAnnotations.put("MedMQ0pc", mq0Median);
		
	}

	private static double rms(int... values) {
		long sqSum = 0;
		long count = 0;
		for (int i = 0; i < values.length; i++) {
			count += values[i];
			sqSum += values[i] * i * i;
		}
		return Math.sqrt(((double) sqSum) / count);
	}

	private static double rms(int[] values, int fromIndex, int toIndex) {
		if (fromIndex < 0 || fromIndex >= values.length)
			throw new IllegalArgumentException(
					"fromIndex must be from 0 to the input length - 1");
		if (toIndex < 0 || toIndex >= values.length)
			throw new IllegalArgumentException(
					"toIndex must be from 0 to the input length");
		if (fromIndex == 0 && toIndex == values.length)
			return rms(values);
		else {
			long sqSum = 0;
			long count = 0;
			for (int i = fromIndex; i < toIndex; i++) {
				count += values[i];
				sqSum += values[i] * i * i;
			}
			return Math.sqrt(((double) sqSum) / count);
		}
	}

	public Set<VCFHeaderLine> getHeaderLines() {
		Set<VCFHeaderLine> result = new HashSet<VCFHeaderLine>();
		result.add(new VCFInfoHeaderLine("DP",1,VCFHeaderLineType.Integer,"total depth at that site across samples"));
		result.add(new VCFInfoHeaderLine("MQ0",1,VCFHeaderLineType.Integer,"total number of reads with mapping quality zero across all samples overlapping this position"));
		result.add(new VCFInfoHeaderLine("MQ", 1, VCFHeaderLineType.Float,"RMS of mapping qualities of reads overlapping this position"));
		result.add(new VCFFormatHeaderLine("DP",1,VCFHeaderLineType.Integer,"total depth at that site across samples"));
		result.add(new VCFFormatHeaderLine("MQ0",1,VCFHeaderLineType.Integer,"total number of reads with mapping quality zero across all samples overlapping this position"));
		result.add(new VCFFormatHeaderLine("MQ", 1, VCFHeaderLineType.Float,"RMS of mapping qualities of reads overlapping this position"));
		result.add(new VCFInfoHeaderLine("MedMQ0pc",1,VCFHeaderLineType.Float,"MQ0 pc median across samples"));
		result.add(new VCFInfoHeaderLine("MedDP",1,VCFHeaderLineType.Float,"DP median across samples"));
		if (walker.normalize) {
			result.add(new VCFInfoHeaderLine("MedCMF",1,VCFHeaderLineType.Float,"Cmf median across samples"));
			result.add(new VCFInfoHeaderLine("Cmf",1,VCFHeaderLineType.Float,"coverage median fold: ratio between the depth at this position and the normalized median"));
			result.add(new VCFInfoHeaderLine("Ccp",1,VCFHeaderLineType.Float,"coverage cumulative probability. e.g. 0.15 means that 15% of positions have lower or equal depth"));
			result.add(new VCFFormatHeaderLine("CMF",1,VCFHeaderLineType.Float,"coverage median fold: ratio between the depth at this position and the normalized median"));
			result.add(new VCFFormatHeaderLine("CCP",1,VCFHeaderLineType.Float,"coverage cumulative probability. e.g. 0.15 means that 15% of positions have lower or equal depth"));
		}
		if (walker.features != null)
			result.add(new VCFInfoHeaderLine("CODING",0,VCFHeaderLineType.Flag,"Indicates that the locus is coding"));
		return result;
	}

}
