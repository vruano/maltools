package net.malariagen.gatk.alignment.test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Formatter;
import java.util.List;

import org.junit.Test;

import net.malariagen.gatk.test.WalkerTest;


public class CompareAlignmentsWalkerTest extends WalkerTest {
	
	@Test
	public void testIndelComparison() throws IOException {
		File reference = new File (getClass().getResource("/testdata/compare/3D7_pm.fa").getFile());
		File features = new File (getClass().getResource("/testdata/compare/3D7_pm.gff").getFile());
		File uniqueness = new File (getClass().getResource("/testdata/compare/3D7_pm.uq").getFile());
		File sampleOne = new File (getClass().getResource("/testdata/compare/PG0060-C.before.bam").getFile());
		File sampleTwo = new File (getClass().getResource("/testdata/compare/PG0060-C.after.bam").getFile());
		File output = File.createTempFile("AGVTestCompAln", ".json");
        Formatter formatter = new Formatter();
        Formatter cmdSpec = formatter.format("-R %s -T CompareAlignments -B:uniqueness,UQN %s -B:features,GFF %s -I %s -I %s -o %s",
        		reference,uniqueness,features,sampleOne,sampleTwo,output);
        // No md5s for now.
        List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(cmdSpec.toString(),md5s);
        executeTest("testCompareAlignments", spec);
	}

}
