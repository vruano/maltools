package net.malariagen.gatk.coverage.test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import java.util.Formatter;

import net.malariagen.gatk.test.WalkerTest;


public class ReadCountsWalkerTest extends WalkerTest {

	@Test
	public void testReadCountsWalker() throws IOException {
		File reference = new File (getClass().getResource("/testdata/reference.fa").getFile());
		File features = new File (getClass().getResource("/testdata/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/sampleTwo.bam").getFile());
//		File coverageDistribution = new File(getClass().getResource("/testdata/sampleOneAndTwo.coverage.json").getFile());
		File output = File.createTempFile("AGVTest", ".vcf");
        Formatter formatter = new Formatter();
//        Formatter cmdSpec = formatter.format("-R %s -A CodingAnnotation -A MinorAlleleCounts -A AbsoluteCounts -T DiploidGenotyper -B:features,GFF %s -I %s -I %s -CvgD %s -o %s",
//        		reference,features,sampleOne,sampleTwo,coverageDistribution,output);
      Formatter cmdSpec = formatter.format("-nt 2 -Q 27 -R %s -A CodingAnnotation -A MinorAlleleCounts -A AbsoluteCounts -T " + "ReadCounts" + " -B:features,GFF %s -I %s -I %s -o %s",
		reference,features,sampleOne,sampleTwo,output);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testReadCountsWalker", spec);
	}
}
