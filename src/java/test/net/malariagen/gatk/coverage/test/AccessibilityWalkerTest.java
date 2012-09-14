package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import java.util.Formatter;

import net.malariagen.gatk.test.WalkerTest;

public class AccessibilityWalkerTest extends WalkerTest {

	@Test
	public void testAccessibilityWalkerNoNormalization() throws IOException {
		File reference = new File (getClass().getResource("/cases/awtest/reference.fa").getFile());
		File features = new File (getClass().getResource("/cases/awtest/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/cases/awtest/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/cases/awtest/sampleTwo.bam").getFile());
		File output = File.createTempFile("cqw", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-R %s -T CoverageQuality -groupBy SM -features %s -I %s -I %s -o %s",
        		reference,features,sampleOne,sampleTwo,output);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        try {
		   executeTest("testAccessibility", spec);
        }
        catch (Throwable t) {
        	t.printStackTrace();
        	fail(t.getMessage());
        }
     }
	
	@Test
	public void testAccessibilityWalkerWithNormalization() throws IOException {
		File reference = new File (getClass().getResource("/cases/awtest/reference.fa").getFile());
		File features = new File (getClass().getResource("/cases/awtest/reference.gff").getFile());
		File sampleOne = new File (getClass().getResource("/cases/awtest/sampleOne.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/cases/awtest/sampleTwo.bam").getFile());
		File cds = new File (getClass().getResource("/cases/awtest/coverage-distribution.json").getFile());
		File output = File.createTempFile("cqw", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-R %s -T CoverageQuality -groupBy SM -features %s -I %s -I %s -o %s -cds %s",
        		reference,features,sampleOne,sampleTwo,output,cds);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        try {
		   executeTest("testAccessibility", spec);
        }
        catch (Throwable t) {
        	t.printStackTrace();
        	fail(t.getMessage());
        }
     }

}
