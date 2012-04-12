package org.broadinstitute.sting.gatk.walkers.genotyper;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


import net.malariagen.gatk.csl.CSLFeature;
import net.malariagen.gatk.genotyper.GenotypeVariantFilterEmitMode;
import net.malariagen.gatk.genotyper.GenotypingContext;
import net.malariagen.gatk.genotyper.MetaArgumentCollection;
import net.malariagen.gatk.genotyper.SnpGenotypingContext;
import net.malariagen.gatk.genotyper.models.GenotypingModel;
import net.malariagen.gatk.genotyper.models.GenotypingModelException;
import net.malariagen.gatk.genotyper.models.GenotypingModelUtils;


import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.SNPGenotypeLikelihoodsCalculationModel;
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
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.MutableGenotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class MetaGenotyperEngine extends UnifiedGenotyperEngine {

	private static final String CANDIATE_SNP_LIST_ROD_NAME = "csl";

	public static final String NO_VARIANT_GT_FILTER = "NoVarGT";
	
	public static final String NO_POLYMORPHIC_GT_FILTER = "NoPolyGT";

	private MetaArgumentCollection metaUAC;

	private PrintWriter baseqDistOut;
	private Logger metaLogger;
	private boolean filterBySnpList = false;

	private VariantAnnotatorEngine metaAnnotationEngine;

	private ThreadLocal<GenotypingModel> genotypingModel = new ThreadLocal<GenotypingModel>();
	
	private int sampleCount;
	
	public MetaGenotyperEngine(GenomeAnalysisEngine toolkit,
			MetaArgumentCollection UAC, Logger logger,
			PrintStream verboseWriter, VariantAnnotatorEngine engine,
			Set<String> samples) {
		super(toolkit, UAC, logger, verboseWriter, engine, samples);
		sampleCount = toolkit.getSamples().size();
		metaUAC = UAC;
		metaLogger = logger;
		metaAnnotationEngine = engine;
		filterBySnpList = checkFilterBySnpList(toolkit);
	}


	private boolean checkFilterBySnpList(GenomeAnalysisEngine toolkit) {
		List<ReferenceOrderedDataSource> rodSources = toolkit.getRodDataSources();
		for (ReferenceOrderedDataSource rod : rodSources) {
			if (rod.getName().equals(CANDIATE_SNP_LIST_ROD_NAME))
				return true;
		}
		return false;
	}
	
	private boolean isCandidatePosition(RefMetaDataTracker rmdt) {
		if (!filterBySnpList) 
			return true;
		
		if (!rmdt.hasROD(CANDIATE_SNP_LIST_ROD_NAME))
			return false;
		List<GATKFeature> features = rmdt.getGATKFeatureMetaData(CANDIATE_SNP_LIST_ROD_NAME, true);
		return features.size() > 0;
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
		genotypingModel.set(result = GenotypingModelUtils.getModelInstance(smodel));
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
	public VariantCallContext calculateLikelihoodsAndGenotypes(
			RefMetaDataTracker tracker, ReferenceContext refContext,
			AlignmentContext rawContext) {
	
		
		if (!isCandidatePosition(tracker)) 
			return null;
		
		
		
		final GenotypingContext gc = buildGenotypingContext(tracker, refContext, rawContext);
		
		PileupElementFilter pueFilter = buildPileupElementFilter(tracker,gc);
		
		rawContext = new AlignmentContext(rawContext.getLocation(), rawContext.getBasePileup().getFilteredPileup(pueFilter),rawContext.getSkippedBases(), rawContext.hasPileupBeenDownsampled());
		if (metaUAC.COVERAGE_AT_WHICH_TO_ABORT > 0
				&& rawContext.size() > metaUAC.COVERAGE_AT_WHICH_TO_ABORT)
			return null;
		
		final GenotypeLikelihoodsCalculationModel.Model model = getCurrentGLModel(
				tracker, refContext, rawContext);
		if (model == null) {
			return (metaUAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES
					&& metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(
					tracker, refContext, null, rawContext) : generateEmptyContext(
							tracker, refContext, null, rawContext));
		}


		Map<String, AlignmentContext> stratifiedContext = getFilteredAndStratifiedContexts(
				metaUAC, refContext, rawContext,
				GenotypeLikelihoodsCalculationModel.Model.SNP);
		if (stratifiedContext == null
				|| gc.getAlleleCount() <= 1) {
			return (metaUAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES
					&& metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(
					tracker, refContext, stratifiedContext, rawContext) : generateEmptyContext(
							tracker, refContext, stratifiedContext, rawContext));
		}

		if (gc.getAlleleCount() <= 1) {
			return (metaUAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES
					&& metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(
					tracker, refContext, stratifiedContext, rawContext) : generateEmptyContext(
							tracker, refContext, stratifiedContext, rawContext));
		}

		GenotypingModel gmodel = getGenotypingModel();
		if (gmodel == null)
			throw new GenotypingModelException(
					"a genotyping model must be specified");
		String[] samples = stratifiedContext.keySet().toArray(new String[stratifiedContext.size()]);
		int sampleCount = samples.length; 
		AlignmentContext[] ac = new AlignmentContext[sampleCount];
		for (int i = 0; i < sampleCount; i++) {
			String s = samples[i];
			ac[i] = stratifiedContext.get(s);
		}

		gmodel.setGenotypingContext(gc);
		Map<String, MutableGenotype> newGenotypes = gmodel
				.calculateGenotypes(stratifiedContext);
		GenomeLoc locus = refContext.getLocus();
		VariantContext newVc = new VariantContext("MG_call", locus.getContig(), locus.getStart(),
				locus.getStop(), gc.getAlleleList(),consolidateGenotypes(newGenotypes).values());

		newVc = metaAnnotationEngine
				.annotateContext(tracker, refContext, stratifiedContext, newVc)
				.iterator().next();
		newVc = filterGenotypeCalls(newVc,gmodel);
		if (newVc == null) return null;
		double negLog10VarQual = gmodel
				.calculateVariantPhredQuality(stratifiedContext,newVc.getGenotypes()) / 10.0;
		if (negLog10VarQual > 99999.99) negLog10VarQual = 99999.99;
		Set<String> filters;
		if (!Double.isNaN(negLog10VarQual)
				&& metaUAC.STANDARD_CONFIDENCE_FOR_CALLING > negLog10VarQual * 10) {
			filters = Collections.singleton(LOW_QUAL_FILTER_NAME);
		} else {
			filters = Collections.emptySet();
		}
		
		newVc = VariantContext.modifyPErrorFiltersAndAttributes(newVc, negLog10VarQual, filters, newVc.getAttributes());
		VariantCallContext result = new VariantCallContext(newVc,
				confidentlyCalled(negLog10VarQual));
		result.setRefBase(gc.getReferenceAllele().getBases()[0]);
		dumpBaseqDistOut(stratifiedContext);
		return result;
	}


	private SnpGenotypingContext buildGenotypingContext(
			RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
		
		byte possibleAlternatives = candidateSnpListAlternative(tracker);
		
		if (possibleAlternatives == 'N')
			return new SnpGenotypingContext(true,refContext,rawContext);
		else
			return new SnpGenotypingContext(true,possibleAlternatives,refContext,rawContext);
	}


	private byte candidateSnpListAlternative(RefMetaDataTracker tracker) {
		// TODO Auto-generated method stub
		List<GATKFeature> tracks = tracker.getGATKFeatureMetaData(CANDIATE_SNP_LIST_ROD_NAME, true);
		for (GATKFeature ft : tracks) {
			Object o = ft.getUnderlyingObject();
			if (o instanceof CSLFeature) {
				CSLFeature c = (CSLFeature)o;
				return c.getAlternative();
			}
		}
		return 'N';
	}


	private PileupElementFilter buildPileupElementFilter(RefMetaDataTracker tracker,
			final GenotypingContext gc) {
		
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
			}};
	}


	private VariantContext filterGenotypeCalls(VariantContext newVc, GenotypingModel gc) {
		Map<String, Genotype> genotypes = newVc.getGenotypes();
		Map<String, Genotype> newGenotypes = null;
		int[] genotypeCounts = new int[gc.getGenotypeCount()];
		int totalDifferentGenotypes = 0;
		for (Map.Entry<String, Genotype> e : genotypes.entrySet()) {
			String s = e.getKey();
			Genotype g = e.getValue();
			
			double gConf = g.hasAttribute("GC") ? g.getAttributeAsDouble("GC") : Double.MIN_VALUE;
			if (g.getNegLog10PError() < metaUAC.MIN_GENOTYPE_QUALITY || gConf < metaUAC.MIN_GENOTYPE_CONFIDENCE) {
				if (newGenotypes == null) 
					newGenotypes = new HashMap<String,Genotype>(genotypes);
				newGenotypes.put(s,Genotype.modifyAlleles(g,Collections.singletonList(Allele.NO_CALL)));
			}
			else {
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
						newVc = addFilterToVariantContext(newVc, NO_POLYMORPHIC_GT_FILTER);
					}
					else if (totalDifferentGenotypes == 0 || genotypeCounts[0] != 0) {
						if (metaUAC.gtVarFilterEmitMode == GenotypeVariantFilterEmitMode.VARIANT)
							return null;
						else 
							newVc = addFilterToVariantContext(newVc, NO_VARIANT_GT_FILTER);
					}
				}
			}
			else if (totalDifferentGenotypes == 0 || genotypeCounts[0] != 0) {
				switch (metaUAC.gtVarFilterEmitMode) {
				case VARIANT:
				case POLYMORPHIC:
					return null;
				default:
					newVc = addFilterToVariantContext(newVc, NO_VARIANT_GT_FILTER);
				}
			}
		}
		return newGenotypes == null ? newVc : VariantContext.modifyGenotypes(newVc, newGenotypes);
	}


	private VariantContext addFilterToVariantContext(VariantContext newVc,
			String newFilter) {
		Set<String> filters = newVc.getFilters();
		if (filters == null || filters.size() == 0) {
			filters = Collections.singleton(newFilter);
		}
		else {
			filters = new HashSet<String>(filters);
			filters.add(newFilter);
		}
		newVc = VariantContext.modifyFilters(newVc, filters);
		return newVc;
	}

	

	private Map<String, Genotype> consolidateGenotypes(
			Map<String, MutableGenotype> newGenotypes) {
		Map<String, Genotype> g = new HashMap<String, Genotype>(
				newGenotypes.size());
		for (Map.Entry<String, MutableGenotype> e : newGenotypes.entrySet()) {
			MutableGenotype mg = e.getValue();
			//if (mg.getNegLog10PError() < metaUAC.MIN_GENOTYPE_QUALITY)
			//	mg.setAlleles(Collections.singletonList(Allele.NO_CALL));
			g.put(e.getKey(), mg.unmodifiableGenotype());
		}
		return g;
	}

	protected boolean confidentlyCalled(double conf) {
		return conf >= metaUAC.STANDARD_CONFIDENCE_FOR_CALLING;
	}



	private GenotypeLikelihoodsCalculationModel.Model getCurrentGLModel(
			final RefMetaDataTracker tracker,
			final ReferenceContext refContext, final AlignmentContext rawContext) {
		if (rawContext.hasExtendedEventPileup()) {
			// todo - remove this code
			if ((metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.INDEL)
					&& (metaUAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES))
				return GenotypeLikelihoodsCalculationModel.Model.INDEL;
		} else {
			// no extended event pileup
			// if we're genotyping given alleles and we have a requested SNP at
			// this position, do SNP
			if (metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) {
				VariantContext vcInput = SNPGenotypeLikelihoodsCalculationModel
						.getSNPVCFromAllelesRod(tracker, refContext, false,
								metaLogger);
				if (vcInput == null)
					return null;

				// todo - no support to genotype MNP's yet
				if (vcInput.isMNP())
					return null;

				if (vcInput.isSNP()) {
					if ((metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP))
						return GenotypeLikelihoodsCalculationModel.Model.SNP;
					else
						// ignore SNP's if user chose INDEL mode
						return null;
				} else if ((vcInput.isIndel() || vcInput.isMixed())
						&& (metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH || metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.INDEL))
					return GenotypeLikelihoodsCalculationModel.Model.INDEL;
			} else {
				// todo - this assumes SNP's take priority when BOTH is
				// selected, should do a smarter way once extended events are
				// removed
				if (metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.BOTH
						|| metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.SNP)
					return GenotypeLikelihoodsCalculationModel.Model.SNP;
				else if (metaUAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.INDEL)
					return GenotypeLikelihoodsCalculationModel.Model.INDEL;
			}
		}
		return null;
	}

	private VariantCallContext generateEmptyContext(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, AlignmentContext> stratifiedContexts,
			AlignmentContext rawContext) {
		VariantContext vc;
		if (metaUAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) {
			VariantContext vcInput = SNPGenotypeLikelihoodsCalculationModel
					.getSNPVCFromAllelesRod(tracker, ref, false, metaLogger);
			if (vcInput == null)
				return null;
			vc = new VariantContext("MG_call", vcInput.getChr(),
					vcInput.getStart(), vcInput.getEnd(), vcInput.getAlleles());
		} else {
			// deal with bad/non-standard reference bases
			if (!Allele.acceptableAlleleBases(new byte[] { ref.getBase() }))
				return null;

			Set<Allele> alleles = new HashSet<Allele>();
			alleles.add(Allele.create(ref.getBase(), true));
			vc = new VariantContext("MG_call", ref.getLocus().getContig(), ref
					.getLocus().getStart(), ref.getLocus().getStart(), alleles);
		}

		if (metaAnnotationEngine != null) {
			// we want to use the *unfiltered* and *unBAQed* context for the
			// annotations
			ReadBackedPileup pileup = null;
			if (rawContext.hasExtendedEventPileup())
				pileup = rawContext.getExtendedEventPileup();
			else if (rawContext.hasBasePileup())
				pileup = rawContext.getBasePileup();
			stratifiedContexts = AlignmentContextUtils
					.splitContextBySampleName(pileup,
							metaUAC.ASSUME_SINGLE_SAMPLE);

			vc = metaAnnotationEngine
					.annotateContext(tracker, ref, stratifiedContexts, vc)
					.iterator().next();
		}

		return new VariantCallContext(vc, ref.getBase(), false);
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
				if (pileup.size() == 0
						&& UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES)
					return null;

				// stratify the AlignmentContext and cut by sample
				stratifiedContexts = AlignmentContextUtils
						.splitContextBySampleName(pileup,
								UAC.ASSUME_SINGLE_SAMPLE);

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
				if (pileup.size() == 0
						&& UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES)
					return null;

				// stratify the AlignmentContext and cut by sample
				stratifiedContexts = AlignmentContextUtils
						.splitContextBySampleName(pileup,
								UAC.ASSUME_SINGLE_SAMPLE);
			}
		} else if (model == GenotypeLikelihoodsCalculationModel.Model.SNP) {

			if (!BaseUtils.isRegularBase(refContext.getBase()))
				return null;

			// stratify the AlignmentContext and cut by sample
			stratifiedContexts = AlignmentContextUtils
					.splitContextBySampleName(rawContext.getBasePileup(),
							UAC.ASSUME_SINGLE_SAMPLE);

			if (!(UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES)) {
				int numDeletions = 0;
				for (final PileupElement p : rawContext.getBasePileup()) {
					if (p.isDeletion()) {
						numDeletions++;
					}
				}
				if (((double) numDeletions)
						/ ((double) rawContext.getBasePileup().size()) > UAC.MAX_DELETION_FRACTION) {
					return null;
				}
			}
		}

		return stratifiedContexts;
	}


	public Collection<? extends VCFHeaderLine> getHeaderLines() {
		return this.getGenotypingModel().getHeaderLines();
	}

}
