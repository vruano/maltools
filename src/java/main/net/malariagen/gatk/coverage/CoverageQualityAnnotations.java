package net.malariagen.gatk.coverage;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.math.IntegerDistributionSet;
import net.malariagen.gatk.math.IntegerSampleDistributionSet;
import net.malariagen.gatk.math.IntegerSequenceDistributionSet;
import net.malariagen.utils.IntRef;

import org.apache.commons.math.stat.descriptive.rank.Median;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;

public class CoverageQualityAnnotations {

	private static final int MAX_MAPPING_QUALITY = 256;

	private CoverageQualityWalker walker;

	private Median medianCalculator = new Median();

	private Map<String, Integer> groupIndex;

	public CoverageQualityAnnotations(CoverageQualityWalker walker) {
		if (walker == null)
			throw new IllegalArgumentException();
		this.walker = walker;
		int nextIdx = 0;
		groupIndex = new LinkedHashMap<String, Integer>(
				walker.groupNames.size());
		for (String gn : walker.groupNames)
			groupIndex.put(gn, nextIdx++);
	}

	public void annotate(RefMetaDataTracker tracker, ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedAlignmentContext,
			Map<String, Object> infoAnnotations,
			Map<String, Map<String, Object>> formatAnnotations) {

		int maxAllMappingQuality;
		int[] allMappingQualityCounts;
		int[][] mappingQualityCounts;
		int[] maxMappingQuality;
		int[] maxForwardStartQuality;
		int[] groupDepth;
		int allDepth;
		double[] groupDoubleAuxiliar;
		double[] groupDoubleAuxiliar2;
		double[] groupDoubleAuxiliar3;

		mappingQualityCounts = new int[groupIndex.size()][MAX_MAPPING_QUALITY + 1];
		maxMappingQuality = new int[groupIndex.size()];
		maxForwardStartQuality = new int[groupIndex.size()];
		maxAllMappingQuality = 0;
		allMappingQualityCounts = new int[MAX_MAPPING_QUALITY + 1];
		groupDepth = new int[groupIndex.size()];
		groupDoubleAuxiliar = new double[groupIndex.size()];
		groupDoubleAuxiliar2 = new double[groupIndex.size()];
		groupDoubleAuxiliar3 = new double[groupIndex.size()];
		allDepth = 0;
		
		IntRef maxAllMappingQuality$ref = new IntRef(maxAllMappingQuality);
		IntRef allDepth$ref = new IntRef(allDepth);
		
		annotate$makeDepthAndMQ0(stratifiedAlignmentContext,
				allMappingQualityCounts, mappingQualityCounts,
				maxMappingQuality, groupDepth, maxForwardStartQuality, maxAllMappingQuality$ref,
				allDepth$ref);
		
		allDepth = allDepth$ref.value;
		maxAllMappingQuality = maxAllMappingQuality$ref.value;

		boolean isCoding = annotate$makeUnormalizeInfoAnnotations(tracker,infoAnnotations,
				maxAllMappingQuality, allMappingQualityCounts, groupDepth,
				allDepth, groupDoubleAuxiliar);
		
		if (walker.normalize) 
			annotate$makeNormalizedAnnotations(ref, infoAnnotations,
					formatAnnotations, allMappingQualityCounts,
					mappingQualityCounts, maxMappingQuality, groupDepth,
					allDepth, groupDoubleAuxiliar, groupDoubleAuxiliar2,
					groupDoubleAuxiliar3, isCoding, maxForwardStartQuality);
		else 
			annotate$makeUnormalizedFormatAnnotations(formatAnnotations,
					mappingQualityCounts, maxMappingQuality, groupDepth,
					groupDoubleAuxiliar,maxForwardStartQuality);
		if (allDepth > 0 && allMappingQualityCounts[0] > 0) {
			double mq0Median = medianCalculator.evaluate(groupDoubleAuxiliar);
			infoAnnotations.put("MedMQ0pc", mq0Median);
		} else {
			infoAnnotations.put("MedMQ0pc", 0.0);
		}
	}

	private void annotate$makeUnormalizedFormatAnnotations(
			Map<String, Map<String, Object>> formatAnnotations,
			int[][] mappingQualityCounts, int[] maxMappingQuality,
			int[] groupDepth, double[] groupDoubleAuxiliar, int[] maxForwardStartQuality) {
		for (String gn : walker.groupNames) {
			int index = groupIndex.get(gn);
			Map<String, Object> annotations = formatAnnotations.get(gn);
			annotations.putAll(annotations);
			annotations.put("MfsQ", maxForwardStartQuality[index]);
			annotations.put("MQ0", mappingQualityCounts[index][0]);
			groupDoubleAuxiliar[index] = groupDoubleAuxiliar[index] == 0 ? 0
					: 100.0 * mappingQualityCounts[index][0]
							/ groupDepth[index];
			annotations.put("DP", groupDepth[index]);
			if (groupDepth[index] != 0)
				annotations.put(
						"MQ",
						rms(mappingQualityCounts[index], 0,
								maxMappingQuality[index] + 1));
		}
	}

	private void annotate$makeNormalizedAnnotations(ReferenceContext ref,
			Map<String, Object> infoAnnotations,
			Map<String, Map<String, Object>> formatAnnotations,
			int[] allMappingQualityCounts, int[][] mappingQualityCounts,
			int[] maxMappingQuality, int[] groupDepth, int allDepth,
			double[] groupDoubleAuxiliar, double[] groupDoubleAuxiliar2,
			double[] groupDoubleAuxiliar3, boolean isCoding, int[] maxForwardStartQuality) {
		IntegerDistribution allCovDist = getAllDistribution(
				walker.coverageDistributionSet, ref, isCoding);
		IntegerDistribution allMq0pcDist = getAllDistribution(
				walker.mq0pcDistributionSet, ref, isCoding);
		if (allCovDist != null) {
			infoAnnotations.put("Cmf", allDepth / allCovDist.median());
			infoAnnotations.put("Ccp",
					allCovDist.cumulativeProbability(allDepth));
		}

		if (allMq0pcDist != null) {
			if (allDepth > 0) {
				infoAnnotations.put("MQ0pcMf", allMappingQualityCounts[0]
						/ (allMq0pcDist.median() * allDepth));
				infoAnnotations.put("MQ0pcCp", allMq0pcDist
						.cumulativeProbability(allMappingQualityCounts[0]
								/ allDepth));
			}
			else {
				infoAnnotations.put("MQ0pcMf", 0);
				infoAnnotations.put("MQ0pcCp", allMq0pcDist
						.cumulativeProbability(0));
				
			}
		}

		for (String gn : walker.groupNames) {
			IntegerDistribution groupCovDist = getSampleDistribution(
					walker.coverageDistributionSet, ref, isCoding, gn);
			IntegerDistribution groupMq0pcDist = getSampleDistribution(
					walker.mq0pcDistributionSet, ref, isCoding, gn);
			int index = groupIndex.get(gn);
			Map<String, Object> annotations = formatAnnotations.get(gn);
			annotations.put("MfsQ", maxForwardStartQuality[index]);
			annotations.put("MQ0", mappingQualityCounts[index][0]);
			groupDoubleAuxiliar[index] = groupDoubleAuxiliar[index] == 0 ? 0
					: 100.0 * mappingQualityCounts[index][0]
							/ groupDepth[index];
			annotations.put("DP", groupDepth[index]);
			if (groupDepth[index] > 0)
				annotations.put(
						"MQ",
						rms(mappingQualityCounts[index], 0,
								maxMappingQuality[index] + 1));
			if (groupCovDist != null) {
				if (groupDepth[index] == 0) {
					annotations.put("CMF", 0);
					annotations.put("CCP",
							groupCovDist.cumulativeProbability(0));
				} else {
					annotations.put("CMF",
							groupDoubleAuxiliar2[index] = groupDepth[index]
									/ groupCovDist.median());
					annotations.put("CCP", groupCovDist
							.cumulativeProbability(groupDepth[index]));
				}
			}
			if (groupMq0pcDist != null) {
				if (groupDepth[index] == 0) {
					annotations.put("MQ0pcMF", 0);
					annotations.put("MQ0pcCP",
							groupMq0pcDist.cumulativeProbability(0));
				} else {
					annotations
							.put("MQ0pcMF",
									groupDoubleAuxiliar3[index] = (100 * mappingQualityCounts[index][0])
											/ (groupMq0pcDist.median() * groupDepth[index]));
					annotations
							.put("MQ0pcCP",
									groupCovDist
											.cumulativeProbability((100 * mappingQualityCounts[index][0])
													/ groupDepth[index]));
				}

			}
		}
		if (allDepth > 0) {
			if (allCovDist != null)
				infoAnnotations.put("MedCMF",
						medianCalculator.evaluate(groupDoubleAuxiliar2));
			if (allMq0pcDist != null && allMappingQualityCounts[0] > 0)
				infoAnnotations.put("MedMQ0pcMF",
						medianCalculator.evaluate(groupDoubleAuxiliar3));
			else if (allMq0pcDist != null)
				infoAnnotations.put("MedMQ0pcMF", 0.0);
		} else {
			if (allCovDist != null)
				infoAnnotations.put("MedCMF", 0.0);
			if (allMq0pcDist != null)
				infoAnnotations.put("MedMQ0pcMF", 0.0);
		}
	}

	private boolean annotate$makeUnormalizeInfoAnnotations(RefMetaDataTracker tracker,
			Map<String, Object> infoAnnotations, int maxAllMappingQuality,
			int[] allMappingQualityCounts, int[] groupDepth, int allDepth,
			double[] groupDoubleAuxiliar) {
		
		infoAnnotations.put("DP", allDepth);
		infoAnnotations.put("MQ0", allMappingQualityCounts[0]);
		if (allDepth > 0)
			infoAnnotations.put("MQ",
					rms(allMappingQualityCounts, 0, maxAllMappingQuality + 1));

		for (int i = 0; i < groupDepth.length; i++)
			groupDoubleAuxiliar[i] = (double) groupDepth[i];

		if (allDepth > 0) {
			double groupMedianDepth = medianCalculator
					.evaluate(groupDoubleAuxiliar);
			infoAnnotations.put("MedDP", groupMedianDepth);
		} else {
			infoAnnotations.put("MedDP", 0.0);
		}

		boolean isCoding = false;
		if (walker.features.isBound()) {
			isCoding = walker.isCoding(tracker);
			if (isCoding)
				infoAnnotations.put("CODING", null);
		}
		if (walker.uniqueness.isBound()) 
			infoAnnotations.put("UQ", walker.uniquenessScore(tracker));
		return isCoding;
	}

	private void annotate$makeDepthAndMQ0(
			Map<String, AlignmentContext> stratifiedAlignmentContext,
			int[] allMappingQualityCounts, int[][] mappingQualityCounts,
			int[] maxMappingQuality, int[] groupDepth, int[] maxForwardStartQuality,
			IntRef maxAllMappingQuality$ref, IntRef allDepth$ref) {
		for (Map.Entry<String, AlignmentContext> gn : stratifiedAlignmentContext
				.entrySet()) {
			int index = groupIndex.get(gn.getKey());
			int[] mappingQualityCount = mappingQualityCounts[index];
			allDepth$ref.value += groupDepth[index] = gn.getValue().getBasePileup()
					.depthOfCoverage();
			for (PileupElement pe : gn.getValue().getBasePileup()) {
				int mq = pe.getMappingQual();
				if (mappingQualityCount[mq]++ == 0
						&& mq > maxMappingQuality[index])
					maxMappingQuality[index] = mq;
				if (allMappingQualityCounts[mq]++ == 0
						&& mq > maxAllMappingQuality$ref.value)
					maxAllMappingQuality$ref.value = mq;
				if (pe.getOffset() == 0 && maxForwardStartQuality[index] < mq)
					maxForwardStartQuality[index] = mq;
			}
		}
	}

	private IntegerDistribution getSampleDistribution(
			IntegerDistributionSet distSet, ReferenceContext ref,
			boolean isCoding, String gn) {
		if (distSet == null)
			return null;
		IntegerSampleDistributionSet gSet = distSet
				.getSampleDistributionSet(gn);
		IntegerSequenceDistributionSet gSeqSet = walker.normalizePerChromosome ? gSet
				.getSequenceDistributionSet(ref.getLocus().getContig()) : gSet
				.getAllSequencesDistributionSet();
		IntegerDistribution gDist = walker.normalizePerCodingStatus ? gSeqSet
				.getCategoryDistribution(isCoding ? LocusCategory.CODING
						: LocusCategory.NON_CODING) : gSeqSet
				.getAllCategoriesDistribution();
		return gDist;
	}

	private IntegerDistribution getAllDistribution(
			IntegerDistributionSet distSet, ReferenceContext ref,
			boolean isCoding) {
		if (distSet == null)
			return null;
		IntegerSampleDistributionSet allSet = distSet
				.getAllSamplesDistributionSet();
		IntegerSequenceDistributionSet allSeqSet = walker.normalizePerChromosome ? allSet
				.getSequenceDistributionSet(ref.getLocus().getContig())
				: allSet.getAllSequencesDistributionSet();
		IntegerDistribution allDist = walker.normalizePerCodingStatus ? allSeqSet
				.getCategoryDistribution(isCoding ? LocusCategory.CODING
						: LocusCategory.NON_CODING) : allSeqSet
				.getAllCategoriesDistribution();

		return allDist;
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

	public static Set<VCFHeaderLine> getHeaderLines(CoverageQualityWalker walker) {
		Set<VCFHeaderLine> result = new HashSet<VCFHeaderLine>();
		result.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer,
				"total depth at that site across samples"));
		result.add(new VCFInfoHeaderLine(
				"MQ0",
				1,
				VCFHeaderLineType.Integer,
				"total number of reads with mapping quality zero across all samples overlapping this position"));
		result.add(new VCFInfoHeaderLine("MQ", 1, VCFHeaderLineType.Float,
				"Mapping Quality RMS of reads overlapping this position"));
		result.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer,
				"Total depth at that site across samples"));
		result.add(new VCFFormatHeaderLine(
				"MQ0",
				1,
				VCFHeaderLineType.Integer,
				"Total number of reads with mapping quality 0 across all samples overlapping this position"));
		result.add(new VCFFormatHeaderLine("MQ", 1, VCFHeaderLineType.Float,
				"RMS of mapping qualities of reads overlapping this position"));
		result.add(new VCFFormatHeaderLine("MfsQ", 1, VCFHeaderLineType.Integer,
				"Maximum Mapping Quality of reads STARTING in the forward strand"));
		result.add(new VCFInfoHeaderLine("MedMQ0pc", 1,
				VCFHeaderLineType.Float, "MQ0 pc median across samples"));
		result.add(new VCFInfoHeaderLine("MedDP", 1, VCFHeaderLineType.Float,
				"Median sample DP"));
		if (!walker.features.isBound())
			result.add(new VCFInfoHeaderLine("CODING",0,VCFHeaderLineType.Flag,"Mark positions that are considered as protein-coding"));
		if (!walker.uniqueness.isBound())
			result.add(new VCFInfoHeaderLine("UQ",1,VCFHeaderLineType.Integer,"Indicates the uniqueness score for that position"));
		
		if (walker.normalize) {
			if (walker.coverageDistributionSet != null) {
				result.add(new VCFInfoHeaderLine("MedCMF", 1,
						VCFHeaderLineType.Float, "Cmf median across samples"));
				result.add(new VCFInfoHeaderLine(
						"Cmf",
						1,
						VCFHeaderLineType.Float,
						"coverage median fold: ratio between the depth at this position and the normalized median"));
				result.add(new VCFInfoHeaderLine(
						"Ccp",
						1,
						VCFHeaderLineType.Float,
						"coverage cumulative probability. e.g. 0.15 means that 15% of positions have lower or equal depth"));
				result.add(new VCFFormatHeaderLine(
						"CMF",
						1,
						VCFHeaderLineType.Float,
						"coverage median fold: ratio between the depth at this position and the normalized median"));
				result.add(new VCFFormatHeaderLine(
						"CCP",
						1,
						VCFHeaderLineType.Float,
						"coverage cumulative probability. e.g. 0.15 means that 15% of positions have lower or equal depth"));
			}
			if (walker.mq0pcDistributionSet != null) {
				result.add(new VCFInfoHeaderLine("MedMQ0pcMF", 1,
						VCFHeaderLineType.Float, "Cmf median across samples"));
				result.add(new VCFInfoHeaderLine(
						"MQ0pcMf",
						1,
						VCFHeaderLineType.Float,
						"Across sample MQ0 percentage median fold: ratio between the MQ0pc at this position and the median"));
				result.add(new VCFInfoHeaderLine(
						"MQ0pcCp",
						1,
						VCFHeaderLineType.Float,
						"Across sample MQ0 percentage cumulative probability. e.g. 0.15 means that 15% of positions have lower or equal MQ0pc"));
				result.add(new VCFFormatHeaderLine(
						"MQ0pcMF",
						1,
						VCFHeaderLineType.Float,
						"MQ0 percentage median fold: ratio between the depth at this position and the median for the sample"));
				result.add(new VCFFormatHeaderLine(
						"MQ0pcCP",
						1,
						VCFHeaderLineType.Float,
						"MQ0 percentage cumulative probability. e.g. 0.15 means that 15% of positions have lower or equal MQ0pc for the sample"));
				result.add(new VCFInfoHeaderLine("MedMQ0pcMF", 1,
						VCFHeaderLineType.Float,
						"MQ0pcMF median across samples at that position"));
			}
		}
		if (walker.features != null)
			result.add(new VCFInfoHeaderLine("CODING", 0,
					VCFHeaderLineType.Flag,
					"Indicates that the locus is coding"));
		return result;
	}

}
