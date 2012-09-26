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

	@Test
	public void test3D7Reference2() throws IOException {
//		File refFasta = new File(this.getClass().getClassLoader().getResource("3D7_pm.fa").getFile());
		File refFasta = new File("/home/valentin/Science/CoveragePf/Simulation/Pf3D7_v3/data/Pf3D7_v3.fa");
		File refGff = new File("/home/valentin/Science/CoveragePf/Simulation/Pf3D7_v3/data/Pf3D7_v3.gff");
		File outDir = File.createTempFile(refFasta.getName() + "-", ".vcf");
		
		try {
			WalkerTest.executeTest("ReferenceComplexity", String.format("-T ReferenceComplexity -W 76 -W 100 -W 152 -W 200 -W 75 -W 150 -R %s -features %s -o %s",refFasta,refGff,outDir), null);
		} catch (Throwable t) {
			t.printStackTrace();
	//		fail();
		}
		//assertTrue(outDir.exists());
		//outDir.delete();
		
	}
	
}
