package net.malariagen.gatk.walker;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@PartitionBy(PartitionType.NONE)
@Requires({ DataSource.READS, DataSource.REFERENCE })
@Downsample(by = DownsampleType.NONE)
public class FragmentDepthAndLengthWalker extends
		ReadWalker<GATKSAMRecord, Region> {

	protected Map<String, Fragment> openFragments = new HashMap<String, Fragment>(
			10000);
	protected PriorityQueue<Fragment> completedFragments = new PriorityQueue<Fragment>();

	protected String currentChromosome;
	protected int lastPosition = 0;
	protected GenomeLocParser locParser;

	@Argument(shortName = "mmq", doc = " minimum mapping quality", required = false)
	public int minimumMappingQuality = 20;

	@Output(shortName = "o", doc = "output file", required = true)
	public File outputFile;

	@Output(shortName = "w", fullName = "where", doc = "file with the positions to output", required = false)
	public RodBinding<Feature> where;

	private Writer output;

	private LocationAwareSeekableRODIterator whereIt;
	private LocationAwareSeekableRODIterator nextIt;

	@Override
	public void initialize() {
		locParser = getToolkit().getGenomeLocParser();
		try {
			output = new FileWriter(outputFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		if (where.isBound()) {
			for (ReferenceOrderedDataSource s : getToolkit()
					.getRodDataSources()) {
				if (s.getName().equals(where.getName())) {
					whereIt = s.seek(null);
					nextIt = s.seek(null);
				}
			}
		}
		super.initialize();
	}

	@Override
	public Region reduceInit() {
		if (locParser == null)
			throw new IllegalStateException();
		return new Region(locParser);
	}

	@Override
	public void onTraversalDone(Region region) {
		try {
			flush(region.loc.getContig(), -1);
			output.flush();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public Region reduce(GATKSAMRecord value, Region sum) {
		try {
			if (value == null)
				return sum;
			if (value.getReadUnmappedFlag())
				return sum;
			if (sum.isMapped()
					&& sum.getContigIndex() != value.getReferenceIndex()) {
	//			flush(sum.loc.getContig(), -1);
	//			currentChromosome = value.getReferenceName();
	//			lastPosition = 0;
				return new Region(locParser, value);
			}
	//		currentChromosome = value.getReferenceName();
			if (completedFragments.size() > 1000)
				flush(currentChromosome, value.getAlignmentStart());
			// lastPosition = value.getAlignmentStart() - 1;
			return sum.merge(value);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private void flush(String chromosome, int to) throws IOException {

		PriorityQueue<Fragment> opens = new PriorityQueue<Fragment>(
				Math.max(openFragments.size(),1), new Comparator<Fragment>() {
					@Override
					public int compare(Fragment arg0, Fragment arg1) {
						return arg0.openningRead.getAlignmentStart()
								- arg1.openningRead.getAlignmentStart();
					}
				});

		for (Fragment f : openFragments.values()) {
			if (!f.openningRead.getReferenceName().equals(chromosome))
				throw new IllegalArgumentException("wrong chromosome " + chromosome + " " + currentChromosome + " " + f.openningRead.getReferenceName());
			opens.add(f);
		}
		if (to <= 0)
			to = getToolkit().getGenomeLocParser().getContigInfo(chromosome)
					.getSequenceLength() + 1;

		Map<Integer, Integer> lengths = new HashMap<Integer, Integer>(
				openFragments.size());
		List<Fragment> toRemove = new LinkedList<Fragment>();
		StringBuffer sb = new StringBuffer("");
		int depth = 0;
		int currentStart = 0;

		for (int i = lastPosition + 1; i < to; i++) {
			boolean changed = i == lastPosition + 1;
			while (!opens.isEmpty()
					&& opens.peek().openningRead.getAlignmentStart() <= i) {
				Fragment f = opens.remove();
				changed = true;
				depth++;
				Integer l = f.openningRead.getInferredInsertSize();
				if (lengths.containsKey(l))
					lengths.put(l, lengths.get(l) + 1);
				else
					lengths.put(l, 1);
			}
			while (!completedFragments.isEmpty()
					&& completedFragments.peek().clossingRead.getAlignmentEnd() < i) {
				Fragment f = completedFragments.remove();
				toRemove.add(f);
				depth--;
				changed = true;
				Integer l = f.openningRead.getInferredInsertSize();
				lengths.put(l, lengths.get(l) - 1);
			}
			if (changed) {
				if (sb.length() > 0) {
					outputBuffer(sb.toString(), chromosome,
							currentStart, i - 1);
					sb.setLength(0);
				}
				currentStart = i;
				createTemplate(chromosome, lengths, sb, depth);
			}
			lastPosition = i;
		}
		if (sb.length() > 0)
			outputBuffer(sb.toString(), chromosome,currentStart,lastPosition);
		for (Fragment f : toRemove) {
			openFragments.remove(f.openningRead.getReadName());
		}
		if (lastPosition == getToolkit().getGenomeLocParser()
				.getContigInfo(chromosome).getSequenceLength()) {
			openFragments.clear();
			completedFragments.clear();
		}
	}

	private void createTemplate(String chromosome,
			Map<Integer, Integer> lengths, StringBuffer sb, int depth) {
		sb.setLength(0);
		sb.append(chromosome).append("\t").append("%d").append("\t")
				.append("%d").append('\t').append(depth).append("\t");
		int depthTest = 0;
		for (Map.Entry<Integer, Integer> e : lengths.entrySet())
			for (int j = 0; j < e.getValue(); j++) {
				sb.append(e.getKey()).append(',');
				depthTest++;
			}
		if (sb.length() > 0)
			sb.setLength(sb.length() - 1);
		if (depthTest != depth)
			throw new IllegalStateException("wrongth depths " + depth
					+ " != " + depthTest);
		sb.append('\n');
	}

	private void outputBuffer(String template, String chr, int start, int end)
			throws IOException {
		if (whereIt == null) {
			output.write(String.format(template, start, end));
			output.flush();
		} else if (!whereIt.hasNext()) {
			// nothing to do.
		} else {
			GenomeLoc nextLoc = whereIt.peekNextLocation();
			GenomeLoc loc = locParser.createGenomeLoc(chr, start, end);
			while (nextLoc != null && nextLoc.isBefore(loc)) {
				if (whereIt.hasNext()) {
					whereIt.next();
					nextLoc = whereIt.peekNextLocation();
				} else {
					nextLoc = null;
				}
			}
			while (nextLoc != null && !nextLoc.isPast(loc)) {
				output.write(String.format(template,
						Math.max(nextLoc.getStart(), start),
						Math.min(nextLoc.getStop(), end)));
				output.flush();
				if (whereIt.hasNext()) {
					whereIt.next();
					nextLoc = whereIt.peekNextLocation();
				} else {
					nextLoc = null;
				}
			}
		}
	}

	@Override
	public GATKSAMRecord map(ReferenceContext ref, GATKSAMRecord read,
			ReadMetaDataTracker metaDataTracker) {

		updateCurrentChromosome(read);
		if (openFragments.containsKey(read.getReadName())) {
			Fragment frag = openFragments.get(read.getReadName());
			frag.clossingRead = read;
			completedFragments.add(frag);
			return read;
		} else {
			if (read.getReadUnmappedFlag())
				return null;
			if (read.getMateUnmappedFlag())
				return null;
			if (read.getMateReferenceIndex() != read.getReferenceIndex())
				return null;
			if (read.getReadNegativeStrandFlag() == read
					.getMateNegativeStrandFlag())
				return null;
			int mateRelOffset = read.getMateAlignmentStart()
					- read.getAlignmentStart();

			if (read.getReadNegativeStrandFlag() && mateRelOffset > 0)
				return null;
			if (!read.getReadNegativeStrandFlag() && mateRelOffset < 0)
				return null;
			if (read.getMappingQuality() < minimumMappingQuality)
				return null;

			if (read.getInferredInsertSize() < 0 || mateRelOffset < 0)
				return null;
			
			if (nextIt != null) {
				GenomeLoc nextLoc = nextIt.hasNext() ? nextIt.peekNextLocation() : null;
				if (nextLoc == null)
					return null;
				GenomeLoc fragLoc = locParser.createGenomeLoc(read.getReferenceName(), read.getAlignmentStart(),read.getAlignmentStart() + read.getInferredInsertSize() - 1);
				while (nextLoc != null && nextLoc.isBefore(fragLoc)) {
					if (nextIt.hasNext()) {
						nextIt.next();
						nextLoc = nextIt.peekNextLocation();
					} else {
						nextLoc = null;
					}
				}
				if (nextLoc == null || fragLoc.isBefore(nextLoc))
					return null;
			}
			
			int mqual = getMateMappingQuality(read);
			if (mqual < minimumMappingQuality)
				return null;
			openFragments.put(read.getReadName(), new Fragment(read));
			return read;
		}
	}

	private void updateCurrentChromosome(GATKSAMRecord read) {
		if (currentChromosome == null)
			currentChromosome = read.getReferenceName();
		else if (!currentChromosome.equals(read.getReferenceName())) {
			try {
				flush(currentChromosome, -1);
				if (!completedFragments.isEmpty() || !openFragments.isEmpty())
					throw new IllegalStateException(""
							+ completedFragments.size() + " "
							+ openFragments.size());
				currentChromosome = read.getReferenceName();
				lastPosition = 0;
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}

	private Map<SAMReaderID, SAMFileReader> samReaders = new HashMap<SAMReaderID, SAMFileReader>();

	private int getMateMappingQuality(GATKSAMRecord read) {
		SAMReaderID srid = getToolkit().getReadsDataSource().getReaderID(read);
		SAMFileReader reader;
		if ((reader = samReaders.get(srid)) == null)
			samReaders.put(srid, reader = new SAMFileReader(getToolkit()
					.getReadsDataSource().getSAMFile(srid)));
		SAMRecordIterator it = reader.query(read.getMateReferenceName(),
				read.getMateAlignmentStart(), read.getMateAlignmentStart() + 1,
				false);
		while (it.hasNext()) {
			SAMRecord candidate = it.next();
			if (candidate.getReadName().equals(read.getReadName())) {
				if (read.getFirstOfPairFlag() == candidate.getFirstOfPairFlag())
					continue;
				it.close();
				return candidate.getMappingQuality();
			}
		}
		it.close();
		throw new IllegalStateException("could not find mate for "
				+ read.getReadName());
	}

}

class Fragment implements Comparable<Fragment> {
	GATKSAMRecord openningRead;
	GATKSAMRecord clossingRead;

	public Fragment(GATKSAMRecord open) {
		openningRead = open;
	}

	public int compareTo(Fragment that) {
		if (this.clossingRead == null && that.clossingRead == null)
			return 0;
		if (this.clossingRead == null && that.clossingRead != null)
			return 1;
		if (this.clossingRead != null && that.clossingRead == null)
			return -1;
		if (this.clossingRead.getReferenceIndex() != that.clossingRead
				.getReferenceIndex())
			return this.clossingRead.getReferenceIndex()
					- that.clossingRead.getReferenceIndex();
		return this.clossingRead.getAlignmentEnd()
				- that.clossingRead.getAlignmentEnd();
	}
}

class Locus implements Comparable<Locus> {

	private GenomeLoc loc;
	Map<String, Map<String, Object>> perGroupData;

	@Override
	public int compareTo(Locus that) {
		return this.loc.compareTo(that.loc);
	}
}

class Region {
	final GenomeLocParser locParser;
	GenomeLoc loc;

	public Region(GenomeLocParser locParser) {
		this.locParser = locParser;
		loc = GenomeLoc.UNMAPPED;
	}

	public boolean isMapped() {
		return !GenomeLoc.isUnmapped(loc);
	}

	public Integer getContigIndex() {
		return loc.getContigIndex();
	}

	public Region(GenomeLocParser locParser, GATKSAMRecord value) {
		this(locParser, locParser.createGenomeLoc(value));
	}

	public Region(GenomeLocParser locParser, GenomeLoc createGenomeLoc) {
		this.locParser = locParser;
		loc = createGenomeLoc;
	}

	public Region merge(GATKSAMRecord value) {
		if (GenomeLoc.isUnmapped(loc))
			return new Region(locParser, value);
		else {
			GenomeLoc valueLoc = locParser.createGenomeLoc(value);
			if (valueLoc.getContigIndex() != loc.getContigIndex())
				throw new IllegalArgumentException(
						"trying to combine different contigs");
			int start = Math.min(loc.getStart(), valueLoc.getStart());
			int stop = Math.max(loc.getStop(), valueLoc.getStop());
			return new Region(locParser, locParser.createGenomeLoc(
					loc.getContig(), start, stop));
		}
	}
}
