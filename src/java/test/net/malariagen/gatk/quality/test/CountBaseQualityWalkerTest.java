package net.malariagen.gatk.quality.test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import java.util.Formatter;

import net.malariagen.gatk.test.WalkerTest;


public class CountBaseQualityWalkerTest extends WalkerTest {

	@Test
	public void testCountBaseQualityWalker() throws IOException {
		File reference = new File (getClass().getResource("/testdata/reference.fa").getFile());
		File features = new File (getClass().getResource("/testdata/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/sampleTwo.bam").getFile());
		File baqOut = File.createTempFile("AGVTestBAQ", ".json");
		File bqOut = File.createTempFile("AGVTestBQ", ".json");
		File mqOut = File.createTempFile("AGVTestMQ", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format(" -R %s -T CountQualities -baqOut %s -bqOut %s -mqOut %s -cs -baq RECALCULATE -B:features,GFF %s -I %s -I %s -log /dev/stderr",
        		reference,baqOut,bqOut,mqOut,features,sampleOne,sampleTwo);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testMappingQualityCounting", spec);
	}

	@Test
	public void testCountBaseQualityWalkerLongData() throws IOException {
		File reference = new File (getClass().getResource("/testdata/reference.fa").getFile());
		File features = new File (getClass().getResource("/testdata/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/AK0009-C.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/AC0001-C.bam").getFile());
		File baqOut = File.createTempFile("AGVTestBAQ", ".json");
		File bqOut = File.createTempFile("AGVTestBQ", ".json");
		File mqOut = File.createTempFile("AGVTestMQ", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format(" -R %s -T CountQualities -baqOut %s -bqOut %s -mqOut %s -cs -baq CALCULATE_AS_NECESSARY -B:features,GFF %s -I %s -I %s -log /dev/stderr",
        		reference,baqOut,bqOut,mqOut,features,sampleOne,sampleTwo);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testMappingQualityCounting", spec);
	}
	
	public static void main(String[] args) throws IOException {
		CountBaseQualityWalkerTest test = new CountBaseQualityWalkerTest();
		//test.testCountMappingQualityWalkerGathererLoad();
		//org.junit.runner.JUnitCore.runClasses(CountCoverageWalkerTest.class);
	}
}
