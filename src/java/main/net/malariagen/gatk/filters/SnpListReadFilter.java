package net.malariagen.gatk.filters;

import java.util.List;

import net.malariagen.gatk.csl.CSLFeature;
import net.malariagen.utils.NucleotideIUPAC;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import scala.actors.threadpool.Arrays;

public class SnpListReadFilter extends ReadFilter {

	@Argument(fullName = "somFilter", shortName = "somFilter", doc = "Indicate whether filter out reads show variation that is not present in the provided candidate list", required = false)
	protected boolean filterReads = false;

	private boolean cslProvided;

	private byte[][] refMask;

	private GenomeAnalysisEngine engine;
	private boolean initialized = false;

	public SnpListReadFilter() {
	}

	public SnpListReadFilter(GenomeAnalysisEngine engine) {
		this();
		initialize(engine);
	}

	@Override
	public void initialize(GenomeAnalysisEngine engine) {
		this.engine = engine;
	}

	private synchronized void reallyIinitialize() {
		if (initialized)
			return;
		initialized = true;
		List<ReferenceOrderedDataSource> sources = engine.getRodDataSources();
		if (sources == null)
			throw new NullPointerException();
		ReferenceDataSource refSource = engine.getReferenceDataSource();
		IndexedFastaSequenceFile reference = refSource.getReference();
		GenomeLocParser glp = new GenomeLocParser(reference);

		List<SAMSequenceRecord> seqRecords = reference.getSequenceDictionary()
				.getSequences();
		refMask = new byte[seqRecords.size()][];
		int nextIdx = 0;
		for (SAMSequenceRecord seqRecord : seqRecords) {
			refMask[nextIdx] = reference.getSequence(
					seqRecord.getSequenceName()).getBases();
			refMask[nextIdx] = Arrays.copyOf(refMask[nextIdx],
					refMask[nextIdx].length);
			nextIdx++;
		}
		Logger logger = Logger.getLogger(SnpListReadFilter.class);
		for (ReferenceOrderedDataSource source : sources) {
			if (source.getName().equals("csl")) {
				if (!cslProvided) {
					cslProvided = true;
					logger.info("CSL: creating reference mask from snp-list file ...");
				}
				parseCslRod(seqRecords, source, glp);
			}
		}
		if (cslProvided)
			logger.info("CSL: finished creating the mask");
	}

	private void parseCslRod(List<SAMSequenceRecord> seqs,
			ReferenceOrderedDataSource source, GenomeLocParser glp) {
		for (int i = 0; i < seqs.size(); i++) {
			SAMSequenceRecord seq = seqs.get(i);
			byte[] mask = refMask[i];
			GenomeLoc gl = glp.createGenomeLoc(seq.getSequenceName(), 1,
					seq.getSequenceLength());
			LocationAwareSeekableRODIterator it = source.seek(gl);
			while (it.hasNext()) {
				RODRecordList rl = it.next();
				for (GATKFeature feature : rl) {
					Object ulo = feature.getUnderlyingObject();
					if (ulo instanceof CSLFeature) {
						CSLFeature csl = (CSLFeature) ulo;
						int pos = csl.getStart();
						mask[pos - 1] = NucleotideIUPAC.fromBases(
								mask[pos - 1], csl.getAlternative())
								.byteValue();
					} else {
						int start = feature.getStart();
						int end = feature.getEnd();
						if (start == end)
							mask[start - 1] = 'N';
						else
							for (int j = start; j <= end; j++)
								mask[j - 1] = 'N';
					}

				}
			}

		}
	}

	/**
	 * Indicate whether any variation in a location would be filtered out
	 * because it appears as invariant in candidate snp list.
	 * 
	 * @param gl
	 *            the target location.
	 * @return {@code false} if the location is not be filtered, {@code true}
	 *         otherwise.
	 */
	public boolean filterOut(GenomeLoc gl) {
		if (!initialized) reallyIinitialize();
		if (!cslProvided)
			return false;

		int refIndex = gl.getContigIndex();
		byte[] mask = refMask[refIndex];
		int start = gl.getStart();
		int stop = gl.getStop();
		if (start == stop)
			return NucleotideIUPAC.isUnambiguous(mask[start - 1]);
		else {
			for (int i = start; i <= stop; i++)
				if (NucleotideIUPAC.isAmbiguous(mask[i - 1]))
					return false;
			return false;
		}
	}

	@Override
	public boolean filterOut(SAMRecord read) {
		if (!initialized) reallyIinitialize();
		if (!cslProvided)
			return false;
		if (read.getReadUnmappedFlag())
			return false;
		int refIndex = read.getReferenceIndex();
		byte[] mask = refMask[refIndex];

		byte[] bases = read.getReadBases();
		for (AlignmentBlock ab : read.getAlignmentBlocks()) {
			int refStart = ab.getReferenceStart() - 1;
			int baseIdx = refStart;
			int readIdx = ab.getReadStart() - 1;
			for (int abLen = ab.getLength(); abLen > 0; abLen--) {
				if (baseIdx >= mask.length) {
					Logger.getLogger(this.engine.getClass()).warn(
							"read alignment block beyond end of chromosome "
									+ read.toString());
					continue;
				}
				if (!NucleotideIUPAC.areCompatible(mask[baseIdx++],
						bases[readIdx++]))
					return true;
			}
		}
		return false;
	}
}
