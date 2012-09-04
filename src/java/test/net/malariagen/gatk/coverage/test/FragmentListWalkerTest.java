package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class FragmentListWalkerTest {

	@Test
	public void test3D7Reference() throws IOException {
	//	File testBam = new File(this.getClass().getClassLoader().getResource("cases/flisttest/data/crosses/bam-files/PG0052-C.bam").getFile());
		File testBam = new File(this.getClass().getClassLoader().getResource("cases/flisttest/PG0051-C.bam").getFile());
		File refFasta = new File(this.getClass().getClassLoader().getResource("cases/flisttest/data/reference/3D7_pm.fa").getFile());
			
		File outDir = File.createTempFile("rctest", ".out.gz");
		outDir.delete();
		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T FragmentList -R %s -I %s -o %s -sort REFERENCE_START -mmq 20 -filter NONE -wf -readGC -mq",
				refFasta,testBam,outDir), null);
		assertTrue(outDir.exists());
//		outDir.delete();
	}
	
}
