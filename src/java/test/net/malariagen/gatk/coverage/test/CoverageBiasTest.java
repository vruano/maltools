package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class CoverageBiasTest {

	@Test
	public void test3D7Reference() throws IOException {
		File testBam = new File(this.getClass().getClassLoader().getResource("cases/cbtest/3D7.bam").getFile());
		File testBam2 = new File(this.getClass().getClassLoader().getResource("cases/cbtest/PG0051-C.bam").getFile());
		File testBam3 = new File(this.getClass().getClassLoader().getResource("cases/cbtest/PG0052-C.bam").getFile());
		File uniqueness = new File(this.getClass().getClassLoader().getResource("cases/cbtest/3D7_pm.uq").getFile());
		File refFasta = new File(this.getClass().getClassLoader().getResource("cases/cbtest/3D7_pm.fa").getFile());
			
		File outDir = File.createTempFile("rctest", ".out");
		outDir.delete();
		outDir.mkdir();
//		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T CoverageBias -R %s -I %s -o %s -fsmmq 10 -groupBy RG -complexity %s ",
//		refFasta,testBam,outDir,complexity),null ); 
		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T CoverageBias -R %s -I %s -L MAL1:1-600000 -groupBy NONE -W 300  -I %s -I %s -o %s -vo %s -mdp 5 -mmq0 0.0  -fsmmq 0  -uqn:uniqueness,UQN %s",
				refFasta,testBam,testBam2,testBam3,outDir,new File(outDir,"site-by-site.vcf"),uniqueness), null);
		assertTrue(outDir.exists());
//		outDir.delete();
	}
	
}
