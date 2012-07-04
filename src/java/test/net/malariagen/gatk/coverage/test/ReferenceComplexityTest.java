package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class ReferenceComplexityTest {

	@Test
	public void test3D7Reference() throws IOException {
		File testBam = new File(this.getClass().getClassLoader().getResource("3D7.bam").getFile());
		File refFasta = new File(this.getClass().getClassLoader().getResource("3D7_pm.fa").getFile());
		File outDir = File.createTempFile("rctest", ".out");
		
		WalkerTest.executeTest("ReferenceComplexity", String.format("-T ReferenceComplexity  -rounding 10 -groupBy WS -R %s -I %s -o %s -fs /home/valentin/Science/CoveragePf/FragmentLengths/gatk-out",refFasta,testBam,outDir), null);
		assertTrue(outDir.exists());
		outDir.delete();
		
	}

}
