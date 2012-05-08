package net.malariagen.gatk.genotyper.models.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Formatter;
import java.util.List;

import net.malariagen.gatk.genotyper.SnpGenotypingContext;
import net.malariagen.gatk.genotyper.models.BernoulliGenotypingModel;
import net.malariagen.gatk.math.Beta;
import net.malariagen.gatk.test.WalkerTest;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.PileupWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.junit.Test;

import scala.actors.threadpool.Arrays;

public class BernoulliGenotypingModelTest extends WalkerTest {

	public static BernoulliGenotypingModel TEST_MODEL_05 = new BernoulliGenotypingModel(
			new SnpGenotypingContext(Allele.create("A", true), Allele.create(
					"T", false)), 0.5);

	@Test
	public void testRB() {
		try {
			double res = Beta.log10(0.9, 18, 200);
			System.err.println(res);
		} catch (RuntimeException e) {
			fail();
		}
		
	}
	
	public void xxxtestX05() {
		assertEquals(2, TEST_MODEL_05.getGenotypeCount());
		List<Allele> zeroGenotype = TEST_MODEL_05.getGenotypeAlleles(0);
		List<Allele> oneGenotype = TEST_MODEL_05.getGenotypeAlleles(1);
		assertEquals(1, zeroGenotype.size());
		assertEquals(zeroGenotype.get(0), Allele.create("A", true));
		assertEquals(1, oneGenotype.size());
		assertEquals(oneGenotype.get(0), Allele.create("T", false));

		List<SAMSequenceRecord> ssr = Collections
				.singletonList(new SAMSequenceRecord("SEQ1", 10));
		SAMSequenceDictionary ssd = new SAMSequenceDictionary(ssr);
		GenomeLocParser glp = new GenomeLocParser(ssd);
		GenomeLoc gl = glp.createGenomeLoc("SEQ1", 1);
		ReferenceContext refContext = new ReferenceContext(glp, gl, (byte) 'A');
		SAMFileHeader sfh = new SAMFileHeader();
		sfh.setSequenceDictionary(ssd);
		SAMRecord sr = new SAMRecord(sfh);
		sr.setAlignmentStart(1);
		List<PileupElement> pel = Arrays.asList(new PileupElement[] {
				new PileupElement(sr, 1), new PileupElement(sr, 2),
				new PileupElement(sr, 3), new PileupElement(sr, 4) });

		SAMFileReader bfr = new SAMFileReader(
				ClassLoader
						.getSystemResourceAsStream("mockups/bernoulli-test/read2.bam"));
		PileupWalker a;
		ReadBackedPileupImpl rbp = new ReadBackedPileupImpl(gl, pel);
		AlignmentContext sample1 = new AlignmentContext(gl, rbp);
	}

	@Test
	public void testWalker() throws IOException {
		File reference = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference.fasta").getFile());
		File sample = new File(getClass().getResource(
				"/mockups/bernoulli-test/read2.bam").getFile());
		File output = File.createTempFile("BernoulliTest", ".vcf");
		Formatter formatter = new Formatter();
		Formatter cmdSpec = formatter
				.format("-R %s -o %s -baseq_do %s.bdo -T MetaGenotyper -I %s -mbq 30 -mmq 0 -out_mode EMIT_ALL_SITES",
						reference, output, output, sample);
		// No md5s for now.
		List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
				cmdSpec.toString(), md5s);
		executeTest("testCoverageCounting", spec);
	}

	@Test
	public void testWalker3() throws IOException {
		File reference = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference.fasta").getFile());
		File sample = new File(getClass().getResource(
				"/mockups/bernoulli-test/read2.bam").getFile());
		File output = File.createTempFile("BernoulliTest", ".vcf");
		Formatter formatter = new Formatter();
		Formatter cmdSpec = formatter
				.format("-R %s -o %s -baseq_do %s.bdo -T MetaGenotyper -I %s -mbq 30 -mmq 0 -out_mode EMIT_ALL_SITES",
						reference, output, output, sample);
		// No md5s for now.
		List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
				cmdSpec.toString(), md5s);
		executeTest("testCoverageCounting", spec);
	}
	
	
	// /home/valentin/src/GV/resources/test/Example3.bam

	@Test
	public void testWalker2() throws IOException {
		File reference = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference/3D7_pm.fa").getFile());
		File uniqness = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference/3D7_pm.uq").getFile());
		File sample = new File(getClass().getResource("/PG0051-C.bam")
				.getFile());
		File sample2 = new File(getClass().getResource("/PG0052-C.bam")
				.getFile());
		File output = File.createTempFile("BernoulliTest", ".vcf");
		Formatter formatter = new Formatter();
		Formatter cmdSpec = formatter
				.format("-R %s -T MetaGenotyper -o %s -baseq_do %s.bdo -I %s -mgc -20 -I %s -out_mode EMIT_ALL_SITES -gem POLYMORPHIC -smodel Bernoulli/0.0001/ -baq RECALCULATE -A AverageBaseQuality -A ReadDepthAndAllelicFractionBySample -A DepthPerAlleleByVariant -A UniquenessScore -A NumberSamplesWithData -B:uniqueness,UQN %s",
						reference, output, output, sample, sample2, uniqness);
		// No md5s for now.
		List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
				cmdSpec.toString(), md5s);
		executeTest("testCoverageCounting", spec);
	}

	@Test
	public void test7G8xGb4() throws IOException {
	//	File baseDir = new File("/data/malariagen/PGV_RD/testbeds/7g8xGb4/som/");
		File baseDir = new File("/data/haldane/malariagen/pfalciparum/crosses/7g8xGb4/som/");
		File _7g8 = new File(baseDir, "PG0084-C.bam");
		File _Gb4 = new File(baseDir, "PG0083-C.bam");

		String[] progenyNames = { "PG0100-C", "PG0101-C", "PG0092-C",
				"PG0112-C", "PG0113-C", "PG0110-C", "PG0111-C", "PG0085-C",
				"PG0095-C", "PG0103-C", "PG0094-C", "PG0087-C", "PG0104-C",
				"PG0107-C", "PG0096-C", "PG0098-C", "PG0106-C", "PG0109-C",
				"PG0088-C", "PG0086-C", "PG0099-C", "PG0093-C", "PG0091-C",
				"PG0105-C", "PG0090-C", "PG0097-C", "PG0108-C", "PG0102-C" };
		File[] bamFiles = new File[progenyNames.length + 2];
		bamFiles[0] = _7g8;
		bamFiles[1] = _Gb4;
		for (int i = 2; i < bamFiles.length; i++)
			bamFiles[i] = new File(baseDir, progenyNames[i - 2] + ".bam");

		StringBuffer sb = new StringBuffer(1000);

		for (File f : bamFiles)
			sb.append("-I ").append(f.toString()).append(' ');

		File reference = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference/3D7_pm.fa").getFile());
		File uniqness = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference/3D7_pm.uq").getFile());
		File output = File.createTempFile("BernoulliTest", ".vcf");
		Formatter formatter = new Formatter();
		Formatter cmdSpec = formatter
				.format("-R %s -T MetaGenotyper -somFilter -o %s -baseq_do %s.bdo -mbq 20 -mmq 20 -mgq 0 -mgc -200 -out_mode EMIT_ALL_SITES -gem EMIT_ALL -smodel Bernoulli/0.01/ -A AverageBaseQuality -A ReadDepthAndAllelicFractionBySample -A DepthPerAlleleByVariant -A UniquenessScore -A NumberSamplesWithData -B:uniqueness,UQN %s %s",
						reference, output, output, uniqness, sb.toString());
		// No md5s for now.
		List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
				cmdSpec.toString(), md5s);
		executeTest("testCoverageCounting", spec);
	}

}
