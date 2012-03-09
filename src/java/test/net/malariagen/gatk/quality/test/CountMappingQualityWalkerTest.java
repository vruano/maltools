package net.malariagen.gatk.quality.test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import java.util.Formatter;

import net.malariagen.gatk.test.WalkerTest;


public class CountMappingQualityWalkerTest extends WalkerTest {

	@Test
	public void testCountMappingQualityWalker() throws IOException {
		File reference = new File (getClass().getResource("/testdata/reference.fa").getFile());
		File features = new File (getClass().getResource("/testdata/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/sampleTwo.bam").getFile());
		File output = File.createTempFile("AGVTestMQ", ".json");
		File output2 = File.createTempFile("AGVTestBSQ", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format(" -R %s -T CountReadQuality -cs -B:features,GFF %s -I %s -I %s -mqOut %s -bqsOut %s",
        		reference,features,sampleOne,sampleTwo,output,output2);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testMappingQualityCounting", spec);
	}

	@Test
	public void testCountMappingQualityWalkerLongData() throws IOException {
		File reference = new File (getClass().getResource("/testdata/reference.fa").getFile());
		File features = new File (getClass().getResource("/testdata/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/AK0009-C.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/AC0001-C.bam").getFile());
		File output = File.createTempFile("AGVTest", ".json");
		File output2 = File.createTempFile("AGVTestBSQ", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-nt 4 -R %s -T CountReadQuality -B:features,GFF %s -I %s -I %s -mqOut %s -bqsOut %s",
        		reference,features,sampleOne,sampleTwo,output,output2);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testCoverageCounting", spec);
	}
	
	public static void main(String[] args) throws IOException {
		CountMappingQualityWalkerTest test = new CountMappingQualityWalkerTest();
		//test.testCountMappingQualityWalkerGathererLoad();
		//org.junit.runner.JUnitCore.runClasses(CountCoverageWalkerTest.class);
	}
}
