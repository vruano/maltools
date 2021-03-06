/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package net.malariagen.gatk.genotyper;

import net.malariagen.gatk.filters.SnpListReadFilter;
import net.sf.picard.filter.NotPrimaryAlignmentFilter;
import net.sf.samtools.SAMReadGroupRecord;

import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.samples.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.MetaGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.*;
import java.io.PrintStream;

/**
 * A variant caller which unifies the approaches of several disparate callers.
 * Works for single-sample and multi-sample data. The user can choose from
 * several different incorporated calculation models.
 */
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
@ReadFilters({ DuplicateReadFilter.class, NotPrimaryAlignmentFilter.class, UnmappedReadFilter.class, SnpListReadFilter.class })
//@Reference(window = @Window(start = -10, stop = 10))
@By(DataSource.REFERENCE)
//@Downsample(by = DownsampleType.BY_SAMPLE, toCoverage = 250)
public class MetaGenotyper extends
		LocusWalker<VariantCallContext, MetaGenotyper.UGStatistics> implements
		TreeReducible<MetaGenotyper.UGStatistics>, AnnotatorCompatibleWalker {

	@ArgumentCollection
	private MetaArgumentCollection UAC = new MetaArgumentCollection();

	// control the output
	@Output(doc = "File to which variants should be written", required = true)
	protected VCFWriter writer = null;

	@Argument(fullName = "debug_file", shortName = "debug_file", doc = "File to print all of the annotated and detailed debugging output", required = false)
	protected PrintStream verboseWriter = null;

	@Argument(fullName = "metrics_file", shortName = "metrics", doc = "File to print any relevant callability metrics output", required = false)
	protected PrintStream metricsWriter = null;

	@Argument(fullName = "annotation", shortName = "A", doc = "One or more specific annotations to apply to variant calls", required = false)
	protected List<String> annotationsToUse = new ArrayList<String>();

	@Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", required = false)
	protected String[] annotationClassesToUse = { "Standard" };

	// the calculation arguments
	private MetaGenotyperEngine MG_engine = null;

	// the annotation engine
	private VariantAnnotatorEngine annotationEngine;

	private List<String> parentSamples;

	// enable deletions in the pileup
	public boolean includeReadsWithDeletionAtLoci() {
		return true;
	}

	// enable extended events for indels
	public boolean generateExtendedEvents() {
		return (UAC.GLmodel != GenotypeLikelihoodsCalculationModel.Model.SNP && UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES);
	}

	/**
	 * Inner class for collecting output statistics from the UG
	 */
	public static class UGStatistics  {
		/** The total number of passes examined -- i.e., the number of map calls */
		long nBasesVisited = 0;

		/**
		 * The number of bases that were potentially callable -- i.e., those not
		 * at excessive coverage or masked with N
		 */
		long nBasesCallable = 0;

		/**
		 * The number of bases called confidently (according to user threshold),
		 * either ref or other
		 */
		long nBasesCalledConfidently = 0;

		/** The number of bases for which calls were emitted */
		long nCallsMade = 0;

		/** The total number of extended events encountered */
		long nExtendedEvents = 0;

		double percentCallableOfAll() {
			return (100.0 * nBasesCallable) / (nBasesVisited - nExtendedEvents);
		}

		double percentCalledOfAll() {
			return (100.0 * nBasesCalledConfidently)
					/ (nBasesVisited - nExtendedEvents);
		}

		double percentCalledOfCallable() {
			return (100.0 * nBasesCalledConfidently) / (nBasesCallable);
		}
	}

	/**
	 * Initialize the samples, output, and genotype calculation model
	 * 
	 **/
	@SuppressWarnings("unchecked")
	public void initialize() {
		// get all of the unique sample names
		// if we're supposed to assume a single sample, do so
		Set<String> samples = new TreeSet<String>();
		getToolkit().getFilters();
		samples = SampleUtils.getSAMFileSamples(getToolkit()
					.getSAMFileHeader());
		getToolkit().getRodDataSources();
		// initialize the verbose writer
		if (verboseWriter != null)
			verboseWriter
					.println("AFINFO\tLOC\tREF\tALT\tMAF\tF\tAFprior\tAFposterior\tNormalizedPosterior");
//		annotationEngine = new VariantAnnotatorEngine(getToolkit(),
//				Arrays.asList(annotationClassesToUse), annotationsToUse);
		annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, Collections.EMPTY_LIST,this,getToolkit());
		MG_engine = new MetaGenotyperEngine(getToolkit(), UAC, logger,
				verboseWriter, annotationEngine, samples);
		if (UAC.parents.size() == 0)
			parentSamples = resolveParentSamplesFromTags(getToolkit());
		else {
			parentSamples = UAC.parents;
			checkParentSamples(parentSamples);
		}
		// initialize the header
		VCFHeader header = new VCFHeader(getHeaderInfo(),samples);
		for (String parent : parentSamples) 
			header.addMetaDataLine(new VCFHeaderLine("PARENT", parent));
		writer.writeHeader(header);
	}
	
	private void checkParentSamples(List<String> parents) {
		Collection<Sample> samples = getToolkit().getSampleDB().getSamples();
		Set<String> parentSoFar = new HashSet<String>(parents.size());
		Set<String> sampleNames = new HashSet<String>(samples.size());
		for (Sample sample : samples)
			sampleNames.add(sample.getID());
		for (String parent : parents) 
			if (!sampleNames.contains(parent)) 
				throw new RuntimeException("Missing parent sample " +  parent);
			else if (parentSoFar.contains(parent))
				throw new RuntimeException("Sample appears more than once in parent list " + parent);
			else 
				parentSoFar.add(parent);
	}

	private List<String> resolveParentSamplesFromTags(GenomeAnalysisEngine toolkit) {

		ArrayList<String> parentSamples = new ArrayList<String>();
		String[] parentPerPosition = new String[2];
		for (SAMReaderID rid : toolkit.getReadsDataSource().getReaderIDs()) {
			List<String> tags = rid.getTags().getPositionalTags();
			for (String tag : tags) {
				if (!tag.toLowerCase().startsWith("parent"))
					continue;
				String parentNum = tag.substring("parent".length()).trim();
				int num = -1;
				if (!parentNum.isEmpty()) {
					try {
					  num = Integer.parseInt(parentNum);
					}
					catch (NumberFormatException e) {
						throw new RuntimeException("not valid parent tag used '" + tag + "'");
					}
					if (num < 0)
						throw new RuntimeException("no valid parent tag used '" + tag + "'");
				}
				if (num >= 0) {
					if (parentPerPosition.length <= num)
						parentPerPosition = Arrays.copyOf(parentPerPosition,num*2);
					for (SAMReadGroupRecord rg : toolkit.getSAMFileHeader(rid).getReadGroups()) {
						if (parentPerPosition[num] != null && !parentPerPosition.equals(rg.getSample()))
							throw new RuntimeException("more than one sample with the same number " + num);
						parentPerPosition[num] = rg.getSample();
					}
				}
				else {
					for (SAMReadGroupRecord rg : toolkit.getSAMFileHeader(rid).getReadGroups()) {
						if (parentSamples.contains(rg.getSample())) 
							continue;
						parentSamples.add(rg.getSample());
					}
				}
			}
		}
		for (int i = 0; i < parentPerPosition.length; i++) {
			String parent = parentPerPosition[i];
			if (parent == null)
				continue;
			if (parentSamples.contains(parent))
				continue;
			int pos = parentSamples.size() < i ? parentSamples.size() : i;
			parentSamples.add(pos,parent);
		}
		return parentSamples;
	}
	
	

	private Set<VCFHeaderLine> getHeaderInfo() {
		Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

		// all annotation fields from VariantAnnotatorEngine
		headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

		// annotation (INFO) fields from UnifiedGenotyper
		if (!UAC.NO_SLOD)
			headerInfo.add(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY,
					1, VCFHeaderLineType.Float, "Strand Bias"));
		headerInfo
				.add(new VCFInfoHeaderLine(VCFConstants.DOWNSAMPLED_KEY, 0,
						VCFHeaderLineType.Flag,
						"Were any of the samples downsampled?"));

		// also, check to see whether comp rods were included
//		List<ReferenceOrderedDataSource> dataSources = getToolkit()
//				.getRodDataSources();
//		for (ReferenceOrderedDataSource source : dataSources) {
//			if (source.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME)) {
//				headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0,
//						VCFHeaderLineType.Flag, "dbSNP Membership"));
//			} else if (source.getName().startsWith(
//					VariantAnnotatorEngine.dbPrefix)) {
//				String name = source.getName().substring(
//						VariantAnnotatorEngine.dbPrefix.length());
//				headerInfo.add(new VCFInfoHeaderLine(name, 0,
//						VCFHeaderLineType.Flag, name + " Membership"));
//			}
//		}

		// FORMAT and INFO fields
		headerInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));

		// FILTER fields
		// if ( UAC.STANDARD_CONFIDENCE_FOR_EMITTING <
		// UAC.STANDARD_CONFIDENCE_FOR_CALLING )
		headerInfo.add(new VCFFilterHeaderLine(
				UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low confidence in the existence of the variant based on read information i.e. qualy of non-reference calls"));

		if (UAC.gtVarFilterEmitMode != GenotypeVariantFilterEmitMode.EMIT_ALL) {
			if (UAC.gtVarFilterEmitMode == GenotypeVariantFilterEmitMode.NO_VARIANT_FILTER) {
				headerInfo
						.add(new VCFFilterHeaderLine(
								MetaGenotyperEngine.NO_VARIANT_GT_FILTER,
								"No variant (non-reference) genotype found with confidence at that site"));
			}
			else if (UAC.gtVarFilterEmitMode == GenotypeVariantFilterEmitMode.NO_POLYMORPHIC_FILTER) {
				headerInfo
						.add(new VCFFilterHeaderLine(
								MetaGenotyperEngine.NO_VARIANT_GT_FILTER,
								"No variant (non-reference) genotype found with confidence at that site; only applies if the number of samples is 1"));
				headerInfo
						.add(new VCFFilterHeaderLine(
								MetaGenotyperEngine.NO_POLYMORPHIC_GT_FILTER,
								"No polymorphic genotypes (multiple non-reference) genotypes found with confidence at that site. Does not apply if number of samples is less than 2"));
			}
		}
		
		headerInfo.addAll(MG_engine.getHeaderLines());

		if (parentSamples.size() > 0) {
			headerInfo.add(new VCFInfoHeaderLine("PGT",-1,VCFHeaderLineType.String,"Parent genotype calls"));
			headerInfo.add(new VCFInfoHeaderLine("SEGREGATING",0,VCFHeaderLineType.Flag,"Indicates that the genotype is segregating considering parent calls"));
		}

		
		return headerInfo;
	}

	/**
	 * Compute at a given locus.
	 * 
	 * @param tracker
	 *            the meta data tracker
	 * @param refContext
	 *            the reference base
	 * @param rawContext
	 *            contextual information around the locus
	 * @return the VariantCallContext object
	 */
	
	public String alleleString(Genotype gt, Iterable<Allele> vcAlleles) {
		if (gt == null || gt.isNoCall())
			return ".";
		List<Allele> alleles = gt.getAlleles();
		if (alleles.size() == 0)
			return ".";
		StringBuffer sb = new StringBuffer(10);
		for (int i = 0; i < alleles.size(); i++) {
			int j = 0;
			Allele a1 = alleles.get(i);
			for (Allele a2 : vcAlleles) {
				if (a2.equals(a1)) {
					sb.append(j).append('/');
					break;
				}
				j++;
			}
		}
		if (sb.length() > 0) sb.setLength(sb.length() - 1);
		return sb.toString();
	}
	
	@Override
	public VariantCallContext map(RefMetaDataTracker tracker,
			ReferenceContext refContext, AlignmentContext rawContext) {
		List<VariantCallContext> vcl = MG_engine.calculateLikelihoodsAndGenotypes(tracker, refContext,
				rawContext);
		if (vcl == null || vcl.size() == 0 || parentSamples.size() == 0) 
			return null;
		
		VariantCallContext vc = vcl.get(0);
		if (vc == null)
			return null;
		StringBuffer sb = new StringBuffer(10);
		String previous = null;
		boolean segregating = false;
		Iterable<Allele> vcAlleles = vc.getAlleles();
		for (String p : parentSamples) {
			Genotype gt = vc.getGenotype(p);
			String curr = alleleString(gt,vcAlleles);
			if (previous != null && !previous.equals(curr) && !curr.equals("."))
				segregating = true;
			if (!curr.equals("."))
				previous = curr;
			sb.append(curr).append(',');
		}
		if (sb.length() > 0)
			sb.setLength(sb.length() -1);
		if (sb.length() == 0)
			return vc;
		Map<String,Object> oldAttrs = vc.getAttributes();
		HashMap<String,Object> attrs = new HashMap<String,Object>(oldAttrs.size()  + 2);
		attrs.putAll(oldAttrs);
		attrs.put("PGT", sb.toString());
		if (segregating)
			attrs.put("SEGREGATING",null);
		VariantContextBuilder vcb = new VariantContextBuilder(vc);
		VariantContext newVc = vcb.attributes(attrs).make();;
		return MG_engine.newVariantCallContext(newVc, (byte)-1,confidentlyCalled(-vc.getLog10PError()));
	}
	
	protected boolean confidentlyCalled(double conf) {
		return conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING;
	}


	public UGStatistics reduceInit() {
		return new UGStatistics();
	}

	public UGStatistics treeReduce(UGStatistics lhs, UGStatistics rhs) {
		lhs.nBasesCallable += rhs.nBasesCallable;
		lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
		lhs.nBasesVisited += rhs.nBasesVisited;
		lhs.nCallsMade += rhs.nCallsMade;
		return lhs;
	}

	public UGStatistics reduce(VariantCallContext value, UGStatistics sum) {
		// we get a point for reaching reduce
		sum.nBasesVisited++;

		// can't call the locus because of no coverage
		if (value == null)
			return sum;

		// A call was attempted -- the base was potentially callable
		sum.nBasesCallable++;

		// the base was confidently callable
		sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;

		// can't make a call here
		if (!value.shouldEmit)
			return sum;

		try {
			// we are actually making a call
			sum.nCallsMade++;
			writer.add(value);
		} catch (IllegalArgumentException e) {
			throw new IllegalArgumentException(
					e.getMessage()
							+ "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
		}

		return sum;
	}

	public void onTraversalDone(UGStatistics sum) {
		logger.info(String.format(
				"Visited bases                                %d",
				sum.nBasesVisited));
		logger.info(String.format(
				"Callable bases                               %d",
				sum.nBasesCallable));
		logger.info(String.format(
				"Confidently called bases                     %d",
				sum.nBasesCalledConfidently));
		logger.info(String.format(
				"%% callable bases of all loci                 %3.3f",
				sum.percentCallableOfAll()));
		logger.info(String.format(
				"%% confidently called bases of all loci       %3.3f",
				sum.percentCalledOfAll()));
		logger.info(String.format(
				"%% confidently called bases of callable loci  %3.3f",
				sum.percentCalledOfCallable()));
		logger.info(String.format(
				"Actual calls made                            %d",
				sum.nCallsMade));

		if (metricsWriter != null) {
			metricsWriter.println(String.format(
					"Visited bases                                %d",
					sum.nBasesVisited));
			metricsWriter.println(String.format(
					"Callable bases                               %d",
					sum.nBasesCallable));
			metricsWriter.println(String.format(
					"Confidently called bases                     %d",
					sum.nBasesCalledConfidently));
			metricsWriter.println(String.format(
					"%% callable bases of all loci                 %3.3f",
					sum.percentCallableOfAll()));
			metricsWriter.println(String.format(
					"%% confidently called bases of all loci       %3.3f",
					sum.percentCalledOfAll()));
			metricsWriter.println(String.format(
					"%% confidently called bases of callable loci  %3.3f",
					sum.percentCalledOfCallable()));
			metricsWriter.println(String.format(
					"Actual calls made                            %d",
					sum.nCallsMade));
		}
	}


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

	@Override
	public RodBinding<VariantContext> getDbsnpRodBinding() {
		return null;
	}
	

}
