package net.malariagen.gatk.math;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import scala.actors.threadpool.Arrays;

public class CombineIntegerDistributionSetTest {

	@Test
	public void testCombineIntegerDistributionSetTest() throws IOException {
		File cov1 = new File(getClass().getResource(
				"/testdata/coverage-distribution/cov-1.json").getFile());
		File cov2 = new File(getClass().getResource(
				"/testdata/coverage-distribution/cov-2.json").getFile());
		File cov3 = new File(getClass().getResource(
				"/testdata/coverage-distribution/cov-3.json").getFile());
		File output = File.createTempFile("AGVTest", ".json");
		String[] args = new String[] { "-I", cov1.toString(), "-I",
				cov2.toString(), "-I", cov3.toString(), "-o", output.toString() };
		System.err.println("Arguments " + Arrays.toString(args));
		int exitCode;
		try {
			exitCode = CombineIntegerDistributionSets.executeWitoutExiting(args);
			assertEquals(exitCode, 0);
		} catch (Exception e) {
			fail("unexepected exception " + e.getMessage());
		}
	}
}
