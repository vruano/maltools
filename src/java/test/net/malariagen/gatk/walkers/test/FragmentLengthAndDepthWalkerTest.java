package net.malariagen.gatk.walkers.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class FragmentLengthAndDepthWalkerTest {
	@Test
	public void test3D7Reference() throws IOException {
		File refFasta = new File("/cluster-home/valentin/Science/pf-crosses/data/genome/sanger/version3/September_2012/Pf3D7_v3.fa");
			
		File outFile = File.createTempFile("fladtest", ".fa");
		outFile.delete();
		File bamFile = new File("/cluster-home/valentin/Science/pf-crosses/data/3d7_v3/bwa_default/alignment_valid/ERR019054.valid.bam");
		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T FragmentDepthAndLength -w /tmp/where.bed -R %s -I %s -o %s ",
				refFasta,bamFile,outFile), null);
		assertTrue(outFile.exists());
	}
}
