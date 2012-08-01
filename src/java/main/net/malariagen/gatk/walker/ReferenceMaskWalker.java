package net.malariagen.gatk.walker;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import net.malariagen.utils.Nucleotide;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils.JexlVCMatchExp;

@PartitionBy(PartitionType.NONE)
@Requires({ DataSource.REFERENCE_BASES, DataSource.REFERENCE_ORDERED_DATA })
public class ReferenceMaskWalker extends RodWalker<ReferenceMaskWalker.Locus, ReferenceMaskWalker.Region> {

	@ArgumentCollection
	protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

	/**
	 * Any variant which overlaps entries from the provided mask rod will be
	 * filtered.
	 */
	@Input(fullName = "mask", doc = "Input ROD mask", required = false)
	public RodBinding<Feature> mask;

	/**
	 * VariantFiltration accepts any number of JEXL expressions (so you can have
	 * two named filters by using --filterName One --filterExpression "X < 1"
	 * --filterName Two --filterExpression "X > 2").
	 */
	@Argument(fullName = "expression", shortName = "E", doc = "One or more expression used with INFO fields to mask", required = false)
	protected ArrayList<String> maskExpressionStrs = new ArrayList<String>();

	@Argument(fullName = "maskNotApplicables", shortName = "maskNAs", doc = "whether we should mask locus where any mask expressions cannot be evaluated", required = false)
	protected boolean maskNotApplicables;

	@Argument(fullName = "maskCharacter", shortName = "maskChr", doc= "charecter to be used for the mask", required = false)
	protected char maskCharacter = 'N';
	
	@Argument(fullName = "lowercaseMask", shortName= "lcMask", doc= "whether a lower case based mask should be applied instead", required = false)
	protected boolean lowercaseMask = false;
	
	@Argument(fullName = "maskNonNucleotides", shortName = "maskNN", doc="whether non-nucleotide bases should be masked as well (e.g. X)", required = false)
	protected boolean maskNonNucleotides = true;
	
	@Argument(fullName = "baseLineWidth", shortName="width", doc="how many bases per line should be output in the output", required = false)
	protected int baseLineWidth = 60;
	
	@Output(shortName = "o", fullName = "fastaOutput", doc = "Fasta mask output", required = true)
	protected File fastaOutputFile;
	
	private List<JexlVCMatchExp> maskExpressions;

	private BufferedWriter fastaOutputWriter;

	static final char MISSING_CHAR = '-';
	
	@Override
	public void initialize() {
		super.initialize();

		ArrayList<String> dummyExpressionNames = new ArrayList<String>(
				maskExpressionStrs.size());
		int nextIdx = 1;
		for (@SuppressWarnings("unused")
		String expr : maskExpressionStrs)
			dummyExpressionNames.add(String.format("mask%d", nextIdx++));
		maskExpressions = VariantContextUtils.initializeMatchExps(
				dummyExpressionNames, maskExpressionStrs);
		if (fastaOutputFile != null)
			try {
				fastaOutputWriter = new BufferedWriter(new FileWriter(fastaOutputFile));
			} catch (IOException e) {
				throw new StingException("could not open fasta out file for writing",e);
			}
	}

	@Override
	public Locus map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		if (tracker == null || ref == null)
			return null;
		boolean maskIt = mustMask(tracker, context);
		Locus result =  new Locus();
		result.loc = ref.getLocus();
		result.maskIt = maskIt;
		result.refNuc = Nucleotide.fromByte(result.refByte = ref.getBase());
		return result;
	}

	private boolean mustMask(RefMetaDataTracker tracker,
			AlignmentContext context) {
	
		
		Collection<VariantContext> vcs;
		try {
			vcs = tracker.getValues(
				variantCollection.variants, context.getLocation());
		}
		catch(RuntimeException e) {
			throw e;
		}
		boolean maskIt = false;
		for (VariantContext vc : vcs) {
			for (VariantContextUtils.JexlVCMatchExp exp : maskExpressions) {
				try {
					if (VariantContextUtils.match(vc, exp))
						maskIt = true;
				} catch (Exception e) {
					// do nothing unless specifically asked to; it just means
					// that the expression isn't defined for this context
					if (maskNotApplicables)
						maskIt = true;
				}
				if (maskIt)
					break;
			}
			if (maskIt)
				break;
		}
		return maskIt;
	}


	@Override
	public Region reduceInit() {
		return new Region();
	}

	@Override
	public Region reduce(Locus value, Region sum) {
		if (value == null)
			return sum;
		if (!sum.accepts(value)) {
			sum.flush(true);
			sum = new Region();
		}
		sum.add(value);
		sum.flush(false);
		return sum;
	}
	
	@Override
	public void onTraversalDone(ReferenceMaskWalker.Region sum) {
		super.onTraversalDone(sum);
		sum.flush(true);
		try {
			fastaOutputWriter.close();
		} catch (IOException e) {
			throw new StingException("problem closing the fasta output");
		}
	}

	class Region {
		GenomeLoc loc;
		StringBuffer sequence;
		int locusCount;
		
		public boolean accepts(Locus locus) {
		
			if (loc == null)
				return true;
			else 
				return loc.onSameContig(locus.loc);
		}
		
		public void add(Locus locus) {
			if (loc == null) {
				loc = locus.loc;
				sequence = new StringBuffer(100);
				sequence.insert(0, missingCharSequence(loc.size()));
			}
			else if (!loc.onSameContig(locus.loc)) 
				throw new StingException("trying to merge region across contings");
			int startDiff = locus.loc.getStart() - loc.getStart();
			if (startDiff < 0) 
				sequence.insert(0,missingCharSequence(-startDiff));
			int endDiff = locus.loc.getStop() - loc.getStop();
			if (endDiff == 1)
				sequence.append(ReferenceMaskWalker.MISSING_CHAR);
			else if (endDiff > 1)
				sequence.append(missingCharSequence(endDiff));
			
			if (startDiff < 0 || endDiff > 0) {
				loc = loc.union(locus.loc);
				startDiff = locus.loc.getStart() - loc.getStart();
			}
			char ch;
			if (locus.maskIt) {
				ch = lowercaseMask ? Character.toLowerCase((char) locus.refByte) : maskCharacter;
			}
			else if (lowercaseMask) {
				ch = Character.toUpperCase((char) locus.refByte);
			}
			else {
				ch = (char) locus.refByte;
			}
			sequence.setCharAt(startDiff,ch);
			locusCount++; 
		}

		private CharSequence missingCharSequence(final int length) {
			return new CharSequence() {

				@Override
				public char charAt(int index) {
					return ReferenceMaskWalker.MISSING_CHAR;
				}

				@Override
				public int length() {
					return length;
				}

				@Override
				public CharSequence subSequence(int start, int end) {
					return missingCharSequence(end - start);
				}
			};
		}
		
		public int size() {
			if (loc == null)
				return 0;
			return loc.size();
		}
		
		public boolean flush(boolean doAll) {
			try {
			if (loc == null)
				return false;
			int i = 0;
			while (i < sequence.length()) {
				String s = sequence.substring(i, Math.min(i + baseLineWidth,sequence.length()));
				if (s.indexOf(MISSING_CHAR) != -1) 
					break;
				if (s.length() < baseLineWidth && !doAll)
					break;
				if (loc.getStart() == 1 && i == 0) {
					fastaOutputWriter.write('>');
					fastaOutputWriter.write(loc.getContig());
					fastaOutputWriter.write('\n');
				}
				fastaOutputWriter.write(s);
				fastaOutputWriter.write('\n');
				i += baseLineWidth;
			}
			int newStart = loc.getStart() + i;
			if (newStart > loc.getStop()) {
				loc = null;
				sequence = null;
			}
			else {
				loc = getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(),newStart,loc.getStop());
				sequence.replace(0, i, "");
				locusCount -= i;
			}
			return i != 0;
			} catch (IOException o) {
				throw new StingException("problems writing the fasta output",o);
			}
		}
	}
	class Locus {
		GenomeLoc loc;
		boolean maskIt;
		Nucleotide refNuc;
		byte refByte;
	}


}


