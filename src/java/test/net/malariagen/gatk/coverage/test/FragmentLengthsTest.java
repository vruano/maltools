package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import net.malariagen.gatk.coverage.FragmentLengths;
import net.malariagen.gatk.coverage.FragmentLengthSummary;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.junit.Test;

public class FragmentLengthsTest {
	
	private static List<String> smIds = Arrays.asList(new String[] { "s1" });
    private static List<GATKSAMReadGroupRecord> rgs = new ArrayList<GATKSAMReadGroupRecord>();
	 
	private static List<String> rgIds = Arrays.asList(new String[] { "rg1" });
	
	private static String seq1Name = "SEQ1";
	private static int seq1Length = 100000;
	private static String seq2Name = "SEQ2";
	private static int seq2Length = 100000;
	
	private static SAMSequenceRecord seq1 = new SAMSequenceRecord(seq1Name, seq1Length);
	private static SAMSequenceRecord seq2 = new SAMSequenceRecord(seq2Name, seq2Length);
	private static List<SAMSequenceRecord> ssr = Arrays.asList(new SAMSequenceRecord[] { seq1, seq2 });
	private static SAMSequenceDictionary ssd = new SAMSequenceDictionary(ssr);
	private static SAMFileHeader sh = new SAMFileHeader();
	
	static {
		try {
		sh.setSequenceDictionary(ssd);
		for (String id : rgIds) {
			GATKSAMReadGroupRecord rg =new GATKSAMReadGroupRecord(id);
			rg.setSample(smIds.get(0));
			sh.addReadGroup(rg);
			rgs.add(rg);
		}
		}
		catch (Throwable t) {
			throw new RuntimeException(t);
		}
	}
	
	private static int readLength = 76;
	
	
	@Test
	public void testCreate() {

		FragmentLengths fl = FragmentLengths.create(smIds, rgIds, 100);
		assertNotNull(fl);
		assert(fl.size() == 0);
	}

	@Test
	public void testMerge() {
		fail("Not yet implemented");
	}

	@Test
	public void testMergeInFragmentLengths() {
		fail("Not yet implemented");
	}

	@Test
	public void testMergeInFragmentLengthArrays() {
		fail("Not yet implemented");
	}

	@Test
	public void testMergeInFragmentLengthFrequencies() {
		fail("Not yet implemented");
	}
	
	@Test
	public void addMappedPair() {
		GATKSAMRecord[] pair = createMappedPair("ReadPair1",seq1Name,10,500);
		FragmentLengths fl = FragmentLengths.create(smIds, rgIds, 100000);
		fl.add(pair[0]);
		assertEquals(0,fl.size());
		fl.add(pair[1]);
		assertEquals(1,fl.size());
		FragmentLengthSummary summary = fl.summary();
		assertNotNull(summary);
		assertNotNull(summary.getSampleFragmentLengths("s1"));
		assertEquals(1,summary.getSampleFragmentLengths("s1").count());
		assertNotNull(summary.getReadGroupFragmentLengths("rg1"));
		assertEquals(1,summary.getReadGroupFragmentLengths("rg1").count());
	}

	public GATKSAMRecord[] createMappedPair(String name, String seq, int start, int end) {
		List<GATKSAMRecord> list = ArtificialSAMUtils.createPair(sh, name, readLength, start, end - readLength + 1, false, false);
		for (GATKSAMRecord r : list)
			r.setReadGroup(rgs.get(0));
		
		return new GATKSAMRecord[] { list.get(0), list.get(1) };
	}
}
