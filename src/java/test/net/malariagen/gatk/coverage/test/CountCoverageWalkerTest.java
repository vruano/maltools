package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.sting.utils.PathUtils;
import org.junit.Test;
import java.util.Formatter;

import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.math.IntegerDistributionSet;
import net.malariagen.gatk.math.IntegerDistributionSetGatherer;
import net.malariagen.gatk.test.WalkerTest;


public class CountCoverageWalkerTest extends WalkerTest {

	@Test
	public void testCoverageDistributionWalker() throws IOException {
		File reference = new File (getClass().getResource("/cases/cctest/reference.fa").getFile());
		File features = new File (getClass().getResource("/cases/cctest/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/cases/cctest/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/cases/cctest/sampleTwo.bam").getFile());
		File output = File.createTempFile("AGVTest", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-R %s -T CoverageDistribution -groupBy SM -features %s -I %s -I %s -o %s",
        		reference,features,sampleOne,sampleTwo,output);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        try {
		   executeTest("testCoverageCounting", spec);
        }
        catch (Throwable t) {
        	t.printStackTrace();
        	fail(t.getMessage());
        }
     }

	@Test
	public void testMQ0PCDistributionWalker() throws IOException {
		File reference = new File (getClass().getResource("/cases/cctest/reference.fa").getFile());
		File features = new File (getClass().getResource("/cases/cctest/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/cases/cctest/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/cases/cctest/sampleTwo.bam").getFile());
		File output = File.createTempFile("mq0pcd", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-R %s -T MQ0PcDistribution -groupBy SM -features %s -I %s -I %s -o %s",
        		reference,features,sampleOne,sampleTwo,output);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        try {
		   executeTest("testCoverageCounting", spec);
        }
        catch (Throwable t) {
        	t.printStackTrace();
        	fail(t.getMessage());
        }
     }
	
	
	@Test
	public void testCountCoverageWalkerLongData() throws IOException {
		File reference = new File (getClass().getResource("/testdata/reference.fa").getFile());
		File features = new File (getClass().getResource("/testdata/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/AK0009-C.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/AC0001-C.bam").getFile());
		File output = File.createTempFile("AGVTest", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-nt 4 -R %s -T CountCoverage -B:features,GFF %s -I %s -I %s -o %s",
        		reference,features,sampleOne,sampleTwo,output);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testCoverageCounting", spec);
	}
	
	@Test
	public void testCountCoverageWalkerGatherer() throws IOException {
		IntegerDistributionSetGatherer gatherer = new IntegerDistributionSetGatherer();
		File inFile = new File(getClass().getResource("/testdata/sampleOneAndTwo.coverage.json").getFile());
		File outFile = File.createTempFile("AGVTest", ".json");
		gatherer.gather(Arrays.asList(inFile,inFile),outFile);
		outFile.deleteOnExit();
		IntegerDistributionSet inCDS = IntegerDistributionSet.read(new FileReader(inFile));
		IntegerDistributionSet outCDS = IntegerDistributionSet.read(new FileReader(outFile));
		IntegerDistribution inAllDS = inCDS.getAllSamplesDistributionSet().getAllSequencesDistributionSet().getAllCategoriesDistribution();
		IntegerDistribution outAllDS = outCDS.getAllSamplesDistributionSet().getAllSequencesDistributionSet().getAllCategoriesDistribution();
		assertEquals(inAllDS.mean(),outAllDS.mean(),0.0001);
		assertEquals(inAllDS.count() * 2 ,outAllDS.count());		
	}
	
	
	@Test
	public void testCountCoverageWalkerGathererLoad() throws IOException {
		File outFile = File.createTempFile("AGVTest", ".json");
		System.err.println(" " + outFile);
		File cov1 = new File(getClass().getResource("/testdata/coverage-distribution/cov-1.json").getFile());
		File covDir = cov1.getParentFile();
		FilenameFilter jsonFilter = new PathUtils.ExtensionFilter("json");
		File[] covFiles = covDir.listFiles(jsonFilter);
		IntegerDistributionSetGatherer gatherer = new IntegerDistributionSetGatherer();		
		gatherer.gather(Arrays.asList(covFiles),outFile);
	}
	public static void main(String[] args) throws IOException {
		CountCoverageWalkerTest test = new CountCoverageWalkerTest();
		test.testCountCoverageWalkerGathererLoad();
	}
}
