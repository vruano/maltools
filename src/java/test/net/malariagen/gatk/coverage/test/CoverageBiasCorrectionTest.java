package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class CoverageBiasCorrectionTest {

	@Test
	public void test3D7Reference() throws IOException {
		File testBam = new File(this.getClass().getClassLoader().getResource("3D7.bam").getFile());
		File testBam2 = new File(this.getClass().getClassLoader().getResource("PG0051-C.bam").getFile());
		File testBam3 = new File(this.getClass().getClassLoader().getResource("PG0052-C.bam").getFile());
		File uniqueness = new File(this.getClass().getClassLoader().getResource("3D7_pm.uq").getFile());
		File refFasta = new File(this.getClass().getClassLoader().getResource("3D7_pm.fa").getFile());
			
		File outDir = File.createTempFile("rctest", ".out");
		outDir.delete();
		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T CoverageBiasCorrection -W 300 -skip 300 -R %s -I %s -o %s -cb /tmp/PG0051-C.notrinuc -fl /home/valentin/Science/CoveragePf/FragmentLengths/gatk-out -groupBy RG",
				refFasta,testBam2,outDir), null);
		assertTrue(outDir.exists());
//		outDir.delete();
	}
	
}
