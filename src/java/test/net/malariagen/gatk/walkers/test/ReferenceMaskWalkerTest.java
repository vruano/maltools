package net.malariagen.gatk.walkers.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import net.malariagen.gatk.test.WalkerTest;

import org.junit.Test;

public class ReferenceMaskWalkerTest {

	@Test
	public void test3D7Reference() throws IOException {
		File testVcf = new File(this.getClass().getClassLoader().getResource("3D7_pm.complexity.vcf").getFile());
		File refFasta = new File(this.getClass().getClassLoader().getResource("3D7_pm.fa").getFile());
			
		String testExpr = "RefMQ < 10.0";
		File outFile = File.createTempFile("rctest", ".fa");
		outFile.delete();
		WalkerTest.executeTest("MappabilityAnalisys", String.format("-T ReferenceMask -R %s -E '%s' -V %s -o %s ",
				refFasta,testExpr,testVcf,outFile), null);
		assertTrue(outFile.exists());
//		outDir.delete();
	}
	
}
