package net.malariagen.gatk.alignment;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.malariagen.gatk.gff.GFFFeature;
import net.malariagen.gatk.uniqueness.UQNFeature;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@ReadFilters({ BadMateFilter.class, UnmappedReadFilter.class,
		NotPrimaryAlignmentFilter.class })
@PartitionBy(PartitionType.CONTIG)
@By(DataSource.READS)
@Requires({ DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES, DataSource.REFERENCE_ORDERED_DATA })
public class CompareAlignmentsWalker extends
		ReadWalker<SAMRecordInfo, AlignmentDifferencesBuffer> {

	@Argument(fullName = "windowSize", shortName = "ws", doc = "Maximum translocation length to report jointly relocation of the same read between alignments (-1 infinite by default)", required = false)
	int windowSize = -1;

	@Argument(fullName = "bufferSize", shortName = "bs", doc = "Maximum number of buffered reads. beyond this point changes are dumped disregarding the window size (default 100000)", required = false)
	int bufferSize = 100000;

	@Output(fullName = "output", shortName = "o", doc = "Where to output the differences", required = false)
	PrintWriter output = new PrintWriter(System.err);
	
	private boolean outputUniqueness = false;
	
	private boolean outputLocusCategory = false;

	private Set<AlignmentDifferencesBuffer> allBuffers;

	private ThreadLocal<AlignmentDifferencesBuffer> bufferByThread;

	SAMReaderID leftReaderID;

	@Override
	public void initialize() {
		super.initialize();
		allBuffers = new HashSet<AlignmentDifferencesBuffer>();
		bufferByThread = new ThreadLocal<AlignmentDifferencesBuffer>();
		Collection<SAMReaderID> readerIds = this.getToolkit()
				.getReadsDataSource().getReaderIDs();
		if (readerIds.size() != 2)
			throw new UserException(
					"for alignment comparisons you must indicate exactly two samples");
		leftReaderID = readerIds.iterator().next();
		readerIDs = new HashMap<SAMFileReader, SAMReaderID>(2);
		List<ReferenceOrderedDataSource> rodSources = getToolkit().getRodDataSources();
		for (ReferenceOrderedDataSource source : rodSources) {
			if (source.getRecordType().equals(GFFFeature.class))
				outputLocusCategory = true;
			else if (source.getRecordType().equals(UQNFeature.class))
				outputUniqueness = true;
		}
	}

	private AlignmentDifferencesBuffer getBuffer() {
		AlignmentDifferencesBuffer buffer = bufferByThread.get();
		if (buffer == null) {
			synchronized (allBuffers) {
				allBuffers.add(buffer = new AlignmentDifferencesBuffer(this));
			}
			bufferByThread.set(buffer);
		}
		return buffer;
	}

	private Map<SAMFileReader, SAMReaderID> readerIDs;

	SAMReaderID getReaderIDFromRecord(SAMRecord r) {
		synchronized (readerIDs) {
			SAMReaderID result = readerIDs.get(r.getFileSource().getReader());
			if (result == null)
				readerIDs.put(
						r.getFileSource().getReader(),
						result = getToolkit().getReadsDataSource().getReaderID(
								r));
			return result;
		}
	}


	@Override
	public AlignmentDifferencesBuffer reduceInit() {
		return getBuffer();
	}

	@Override
	public AlignmentDifferencesBuffer reduce(SAMRecordInfo value,
			AlignmentDifferencesBuffer sum) {
		sum.add(value);
		return sum;
	}
	
	@Override
	public void onTraversalDone(AlignmentDifferencesBuffer sum) {
		super.onTraversalDone(sum);
		sum.flush();
	}

	@Override
	public SAMRecordInfo map(ReferenceContext ref, GATKSAMRecord read,
			ReadMetaDataTracker metaDataTracker) {
		if (ref == null)
			throw new IllegalArgumentException("");
		return new SAMRecordInfo(read,ref,metaDataTracker,outputUniqueness,outputLocusCategory);
	}
	
	

}
