package org.broadinstitute.sting.gatk.walkers.genotyper;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.genotyper.GenotypeVariantFilterEmitMode;
import net.malariagen.gatk.genotyper.GenotypingContext;
import net.malariagen.gatk.genotyper.MetaArgumentCollection;
import net.malariagen.gatk.genotyper.SnpGenotypingContext;
import net.malariagen.gatk.genotyper.models.GenotypingModel;
import net.malariagen.gatk.genotyper.models.GenotypingModelException;
import net.malariagen.gatk.genotyper.models.GenotypingModelUtils;
import net.malariagen.utils.NucleotideIUPAC;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElementFilter;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

public class MetaGenotyperEngine extends UnifiedGenotyperEngine {

	private static final String NO_ID = ".";

	public class MyVariantContext extends VariantContext {

		public MyVariantContext(String source, String contig, long start,
				long stop, Collection<Allele> alleles,
				GenotypesContext genotypes, double log10PError,
				Set<String> filters, Map<String, Object> attributes) {
			super(source, NO_ID, contig, start, stop, alleles, genotypes,
					log10PError, filters, attributes, NucleotideIUPAC.GAP
							.byteValue(), EnumSet.noneOf(Validation.class));
		}

		public MyVariantContext(VariantContext other) {
			super(other);
		}

		@Override
		public boolean isSNP() {
			return true;
		}

		@Override
		public VariantContext.Type getType() {
			return VariantContext.Type.SNP;
		}

		@Override
		public Allele getAltAlleleWithHighestAlleleCount() {
			if (this.alleles.size() >= 2)
				return super.getAltAlleleWithHighestAlleleCount();
			else {
				Allele ref = this.getReference();
				switch (NucleotideIUPAC.fromBase(ref.getBases()[0])) {
				case A:
					return Allele.create((byte) 'T', false);
				case T:
					return Allele.create((byte) 'A', false);
				case G:
					return Allele.create((byte) 'C', false);
				case C:
					return Allele.create((byte) 'G', false);
				}
				throw new IllegalStateException("invalid reference allele");
			}
		}

	}

	// private static final String CANDIATE_SNP_LIST_ROD_NAME = "csl";

	public static final String NO_VARIANT_GT_FILTER = "NoVarGT";

	public static final String NO_POLYMORPHIC_GT_FILTER = "NoPolyGT";

	private MetaArgumentCollection metaUAC;

	private PrintWriter baseqDistOut;
	// private Logger metaLogger;
	// private boolean filterBySnpList = false;

	private VariantAnnotatorEngine metaAnnotationEngine;

	private ThreadLocal<GenotypingModel> genotypingModel = new ThreadLocal<GenotypingModel>();

	private int sampleCount;
	
	private Set<String> sampleNames;

	public MetaGenotyperEngine(GenomeAnalysisEngine toolkit,
			MetaArgumentCollection UAC, Logger logger,
			PrintStream verboseWriter, VariantAnnotatorEngine engine,
			Set<String> samples) {
		super(toolkit, UAC, logger, verboseWriter, engine, samples, 1);
		sampleCount = samples.size();
		sampleNames = samples;
		metaUAC = UAC;
		// metaLogger = logger;
		metaAnnotationEngine = engine;
		// filterBySnpList = checkFilterBySnpList(toolkit);
	}

	private boolean isCandidatePosition(RefMetaDataTracker rmdt) {
		return true;
	}

	public MetaGenotyperEngine(GenomeAnalysisEngine toolkit,
			MetaArgumentCollection UAC) {
		super(toolkit, UAC);
		metaUAC = UAC;
	}

	public GenotypingModel getGenotypingModel() {
		GenotypingModel result = genotypingModel.get();
		if (result != null)
			return result;
		String smodel = metaUAC.smodel;
		if (smodel == null)
			return null;
		genotypingModel.set(result = GenotypingModelUtils
				.getModelInstance(smodel));
		return result;
	}

	public PrintWriter getBaseqDistOutWriter() {
		if (baseqDistOut != null)
			return baseqDistOut;
		if (metaUAC.baseqDistOut == null)
			return null;
		try {
			return baseqDistOut = new PrintWriter(new FileWriter(
					metaUAC.baseqDistOut));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Compute full calls at a given locus. Entry point for engine calls from
	 * the UnifiedGenotyper.
	 * 
	 * @param tracker
	 *            the meta data tracker
	 * @param refContext
	 *            the reference base
	 * @param rawContext
	 *            contextual information around the locus
	 * @return the VariantCallContext object
	 */
	@Override
	public List<VariantCallContext> calculateLikelihoodsAndGenotypes(
			RefMetaDataTracker tracker, ReferenceContext refContext,
			AlignmentContext rawContext) {

		if (!isCandidatePosition(tracker))
			return null;

		GenotypingContext gc = buildGenotypingContext(tracker, refContext,
				rawContext);

		PileupElementFilter pueFilter = buildPileupElementFilter(tracker, gc);
		rawContext = new AlignmentContext(rawContext.getLocation(), rawContext
				.getBasePileup().getFilteredPileup(pueFilter),
				rawContext.getSkippedBases(),
				rawContext.hasPileupBeenDownsampled());
		gc = buildGenotypingContext(tracker, refContext, rawContext);

		// TODO removed during migration to new GATK as variable has disapeared.
		// if (metaUAC.COVERAGE_AT_WHICH_TO_ABORT > 0
		// && rawContext.size() > metaUAC.COVERAGE_AT_WHICH_TO_ABORT)
		// return null;

		// final GenotypeLikelihoodsCalculationModel.Model model =
		// getCurrentGLModel(
		// tracker, refContext, rawContext);
		// if (model == null) {
		// return Collections
		// .singletonList((metaUAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES
		// && metaUAC.GenotypingMode ==
		// GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
		// ? generateEmptyContext(
		// tracker, refContext, null, rawContext) : null));
		// }

		Map<String, AlignmentContext> stratifiedContext = getFilteredAndStratifiedContexts(
				metaUAC, refContext, rawContext,
				GenotypeLikelihoodsCalculationModel.Model.SNP);
		if (stratifiedContext == null || rawContext.size() == 0) {
			return Collections
					.singletonList((metaUAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES
							&& metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(
							tracker, refContext, stratifiedContext, rawContext)
							: null));
		}

		GenotypingModel gmodel = getGenotypingModel();
		if (gmodel == null)
			throw new GenotypingModelException(
					"a genotyping model must be specified");
		String[] samples = stratifiedContext.keySet().toArray(
				new String[stratifiedContext.size()]);
		int sampleCount = samples.length;
		AlignmentContext[] ac = new AlignmentContext[sampleCount];
		for (int i = 0; i < sampleCount; i++) {
			String s = samples[i];
			ac[i] = stratifiedContext.get(s);
		}

		gmodel.setGenotypingContext(gc);
		GenotypesContext newGenotypes = gmodel
				.calculateGenotypes(stratifiedContext);
		GenomeLoc locus = refContext.getLocus();
		VariantContext newVc = new MyVariantContext("MG_call",
				locus.getContig(), locus.getStart(), locus.getStop(),
				gc.getAlleleList(), newGenotypes, -1, null, null);

		newVc = metaAnnotationEngine.annotateContext(tracker, refContext,
				stratifiedContext, newVc);

		newVc = filterGenotypeCalls(newVc, gmodel);
		if (newVc == null)
			return null;
		double negLog10VarQual = gmodel.calculateVariantPhredQuality(
				stratifiedContext, newVc.getGenotypes()) / 10.0;
		if (negLog10VarQual > 99999.99)
			negLog10VarQual = 99999.99;
		Set<String> filters;
		if (!Double.isNaN(negLog10VarQual)
				&& metaUAC.STANDARD_CONFIDENCE_FOR_EMITTING > negLog10VarQual * 10)
			return null;
		if (!Double.isNaN(negLog10VarQual)
				&& metaUAC.STANDARD_CONFIDENCE_FOR_CALLING > negLog10VarQual * 10)
			filters = Collections.singleton(LOW_QUAL_FILTER_NAME);
		else
			filters = Collections.emptySet();

		newVc = crossAnnotations(newVc);
		VariantContextBuilder vcb = new VariantContextBuilder(newVc);
		vcb.log10PError(-negLog10VarQual);
		vcb.filters(filters);
		newVc = vcb.make();
		VariantCallContext result = new VariantCallContext(newVc,
				confidentlyCalled(negLog10VarQual));
		dumpBaseqDistOut(stratifiedContext);
		return Collections.singletonList(result);
	}

	private VariantContext crossAnnotations(VariantContext vc) {
		return vc;
	}

	private SnpGenotypingContext buildGenotypingContext(
			RefMetaDataTracker tracker, ReferenceContext refContext,
			AlignmentContext rawContext) {

		byte possibleAlternatives = candidateSnpListAlternative(tracker);

		if (possibleAlternatives == 'N')
			return new SnpGenotypingContext(true, refContext, rawContext);
		else
			return new SnpGenotypingContext(true, possibleAlternatives,
					refContext, rawContext);
	}

	private byte candidateSnpListAlternative(RefMetaDataTracker tracker) {
		// List<GATKFeature> tracks =
		// tracker.getGATKFeatureMetaData(CANDIATE_SNP_LIST_ROD_NAME, true);
		// for (GATKFeature ft : tracks) {
		// Object o = ft.getUnderlyingObject();
		// if (o instanceof CSLFeature) {
		// CSLFeature c = (CSLFeature)o;
		// return c.getAlternative();
		// }
		// }
		return 'N';
	}

	private PileupElementFilter buildPileupElementFilter(
			RefMetaDataTracker tracker, final GenotypingContext gc) {

		return new PileupElementFilter() {
			@Override
			public boolean allow(PileupElement pileupElement) {
				byte b = pileupElement.getBase();
				if (gc.getAlleleIndex(b) < 0)
					return false;
				if (pileupElement.getMappingQual() < metaUAC.MIN_MAPPING_QUALTY_SCORE)
					return false;
				if (pileupElement.getQual() < metaUAC.MIN_BASE_QUALTY_SCORE)
					return false;
				return true;
			}
		};
	}

	private VariantContext filterGenotypeCalls(VariantContext newVc,
			GenotypingModel gc) {
		GenotypesContext genotypes = newVc.getGenotypes();
		int[] genotypeCounts = new int[gc.getGenotypeCount()];
		boolean modified = false;
		int totalDifferentGenotypes = 0;
		for (int i = 0; i < genotypes.size(); i++) {
			Genotype g = genotypes.get(i);
			double gConf = g.hasAttribute("GC") ? g.getAttributeAsDouble("GC",
					0) : Double.MIN_VALUE;
			if (-g.getLog10PError() < metaUAC.MIN_GENOTYPE_QUALITY
					|| gConf < metaUAC.MIN_GENOTYPE_CONFIDENCE) {
				modified = true;
				genotypes.set(
						i,
						Genotype.modifyAlleles(g,
								Collections.singletonList(Allele.NO_CALL)));
			} else {
				int gIdx = gc.getGenotypeIndex(g.getAlleles());
				if (gIdx != -1 && genotypeCounts[gIdx]++ == 0)
					totalDifferentGenotypes++;
			}
		}
		if (metaUAC.gtVarFilterEmitMode != GenotypeVariantFilterEmitMode.EMIT_ALL) {
			if (sampleCount > 1) {
				if (totalDifferentGenotypes < 2) {
					if (metaUAC.gtVarFilterEmitMode == GenotypeVariantFilterEmitMode.POLYMORPHIC)
						return null;
					else if (metaUAC.gtVarFilterEmitMode == GenotypeVariantFilterEmitMode.NO_POLYMORPHIC_FILTER) {
						newVc = addFilterToVariantContext(newVc,
								NO_POLYMORPHIC_GT_FILTER);
					} else if (totalDifferentGenotypes == 0
							|| genotypeCounts[0] != 0) {
						if (metaUAC.gtVarFilterEmitMode == GenotypeVariantFilterEmitMode.VARIANT)
							return null;
						else
							newVc = addFilterToVariantContext(newVc,
									NO_VARIANT_GT_FILTER);
					}
				}
			} else if (totalDifferentGenotypes == 0 || genotypeCounts[0] != 0) {
				switch (metaUAC.gtVarFilterEmitMode) {
				case VARIANT:
				case POLYMORPHIC:
					return null;
				default:
					newVc = addFilterToVariantContext(newVc,
							NO_VARIANT_GT_FILTER);
				}
			}
		}
		if (!modified)
			return newVc;
		VariantContextBuilder vcb = new VariantContextBuilder(newVc);
		vcb.genotypes(genotypes);
		return vcb.make();
	}

	private VariantContext addFilterToVariantContext(VariantContext newVc,
			String newFilter) {
		Set<String> filters = newVc.getFilters();
		if (filters == null || filters.size() == 0) {
			filters = Collections.singleton(newFilter);
		} else {
			filters = new HashSet<String>(filters);
			filters.add(newFilter);
		}
		VariantContextBuilder vcb = new VariantContextBuilder(newVc);
		vcb.filters(filters);
		return vcb.make();
	}

	protected boolean confidentlyCalled(double conf) {
		return conf >= metaUAC.STANDARD_CONFIDENCE_FOR_CALLING;
	}

	private VariantCallContext generateEmptyContext(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts,
			AlignmentContext rawContext) {
		VariantContext vc;

		// TODO need to handle this possibility?.
		if (metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES)
			throw new UnsupportedOperationException(
					"do not support genotyping on provided alleles");
		if (!Allele.acceptableAlleleBases(new byte[] { ref.getBase() }))
			return null;
		Set<Allele> alleles = new HashSet<Allele>();
		alleles.add(Allele.create(ref.getBase(), true));
		VariantContextBuilder vcb = new VariantContextBuilder();
		GenomeLoc loc = ref.getLocus();
		vcb.alleles(alleles)
				.loc(loc.getContig(), loc.getStart(), loc.getStop())
				.source("MG_call");
		vc = vcb.make();

		if (metaAnnotationEngine != null) {
			ReadBackedPileup pileup = null;
			if (rawContext.hasExtendedEventPileup())
				pileup = rawContext.getExtendedEventPileup();
			else if (rawContext.hasBasePileup())
				pileup = rawContext.getBasePileup();
			stratifiedContexts = addMissingSamples(ref.getLocus(),AlignmentContextUtils
					.splitContextBySampleName(pileup));
			vc = metaAnnotationEngine.annotateContext(tracker, ref,
					stratifiedContexts, vc);
		}
		return new VariantCallContext(vc, false);
	}

	private Map<String, AlignmentContext> addMissingSamples(GenomeLoc loc,
			Map<String, AlignmentContext> splitContextBySampleName) {
		for (String s : sampleNames) {
			if (splitContextBySampleName.containsKey(s)) continue;
			AlignmentContext ac = new AlignmentContext(loc, new ReadBackedPileupImpl(loc));
			splitContextBySampleName.put(s,ac);
		}
		return splitContextBySampleName;
	}

	private void dumpBaseqDistOut(
			Map<String, AlignmentContext> stratifiedContext) {
		PrintWriter pw = getBaseqDistOutWriter();
		if (pw == null)
			return;
		StringBuffer sb = new StringBuffer();
		for (AlignmentContext ctx : stratifiedContext.values()) {
			ReadBackedPileup rbp = ctx.getBasePileup();
			for (PileupElement pe : rbp) {
				int x = pe.getQual();
				sb.append(x).append(',');
			}
			if (sb.length() > 0)
				sb.setLength(sb.length() - 1);
		}
		pw.println(sb.toString());
		pw.flush();
	}

	private Map<String, AlignmentContext> getFilteredAndStratifiedContexts(
			UnifiedArgumentCollection UAC, ReferenceContext refContext,
			AlignmentContext rawContext,
			final GenotypeLikelihoodsCalculationModel.Model model) {

		Map<String, AlignmentContext> stratifiedContexts = null;

		if (model == GenotypeLikelihoodsCalculationModel.Model.INDEL) {

			if (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) {
				// regular pileup in this case
				ReadBackedPileup pileup = rawContext.getBasePileup();

				// don't call when there is no coverage
				if (pileup.isEmpty()
						&& UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES)
					return null;

				// stratify the AlignmentContext and cut by sample
				stratifiedContexts = addMissingSamples(refContext.getLocus(),AlignmentContextUtils
						.splitContextBySampleName(pileup));

			} else {
				// todo - tmp will get rid of extended events so this wont be
				// needed
				if (!rawContext.hasExtendedEventPileup())
					return null;
				ReadBackedExtendedEventPileup rawPileup = rawContext
						.getExtendedEventPileup();

				// filter the context based on min mapping quality
				ReadBackedExtendedEventPileup pileup = rawPileup;

				// don't call when there is no coverage
				if (pileup.isEmpty()
						&& UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES)
					return null;

				// stratify the AlignmentContext and cut by sample
				stratifiedContexts = addMissingSamples(refContext.getLocus(),AlignmentContextUtils
						.splitContextBySampleName(pileup));
			}
		} else if (model == GenotypeLikelihoodsCalculationModel.Model.SNP) {

			if (!BaseUtils.isRegularBase(refContext.getBase()))
				return null;

			// stratify the AlignmentContext and cut by sample
			stratifiedContexts = addMissingSamples(refContext.getLocus(),AlignmentContextUtils
					.splitContextBySampleName(rawContext.getBasePileup()));

			if (!(UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES)) {
				int numDeletions = 0;
				for (final PileupElement p : rawContext.getBasePileup()) {
					if (p.isDeletion()) {
						numDeletions++;
					}
				}
				if (((double) numDeletions)
						/ ((double) rawContext.getBasePileup()
								.depthOfCoverage()) > UAC.MAX_DELETION_FRACTION) {
					return null;
				}
			}
		}

		return stratifiedContexts;
	}

	public Collection<? extends VCFHeaderLine> getHeaderLines() {
		return this.getGenotypingModel().getHeaderLines();
	}

	public VariantCallContext newVariantCallContext(VariantContext newVc,
			byte ref, boolean confidentlyCalled) {
		VariantCallContext result = new VariantCallContext(newVc,
				confidentlyCalled);
		return result;
	}

}
