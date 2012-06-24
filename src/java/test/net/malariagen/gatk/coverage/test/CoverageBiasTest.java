package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class CoverageBiasTest {

	@Test
	public void test3D7Reference() throws IOException {
		File testBam = new File(this.getClass().getClassLoader().getResource("3D7.bam").getFile());
//		File testBam2 = new File(this.getClass().getClassLoader().getResource("PG0051-C.bam").getFile());
//		File testBam3 = new File(this.getClass().getClassLoader().getResource("PG0052-C.bam").getFile());
//		File uniqueness = new File(this.getClass().getClassLoader().getResource("3D7_pm.uq").getFile());
		File refFasta = new File(this.getClass().getClassLoader().getResource("3D7_pm.fa").getFile());
		File complexity = new File(this.getClass().getClassLoader().getResource("3D7_pm-complexity.vcf").getFile());
			
		File outDir = File.createTempFile("rctest", ".out");
		outDir.delete();
		outDir.mkdir();
		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T CoverageBias -R %s -I %s -o %s -fsmmq 10 -groupBy RG -complexity %s ",
		refFasta,testBam,outDir,complexity),null ); 
//		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T CoverageBias -R %s -I %s -I %s -o %s -fsmmq 10 -groupBy RG -uqn:uniqueness,UQN %s ",
//				refFasta,testBam2,testBam3,outDir,uniqueness), null);
		assertTrue(outDir.exists());
		outDir.delete();
	}
	
}
