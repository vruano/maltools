package net.malariagen.gatk.coverage.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.Test;

import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.gatk.math.IntegerDistributionSet;

public class CountDistributionTest {

	@Test
	public void testDistributionRead() throws IOException {
		FileReader r = new FileReader(getClass().getResource("/testdata/sampleOneAndTwo.coverage.json").getFile());
		IntegerDistributionSet dset = IntegerDistributionSet.read(r);
		assertNotNull(dset);
		r.close();	
	}
	
	@Test
	public void testDistributionWrite() throws IOException {
		File orig = new File(getClass().getResource("/testdata/sampleOneAndTwo.coverage.json").getFile());
		FileReader r = new FileReader(orig);
		IntegerDistributionSet dset = IntegerDistributionSet.read(r);
		File temp = File.createTempFile("AGVTest",".json");
		System.err.println(temp);
		FileWriter w = new FileWriter(temp);
		dset.write(w);
		orig.compareTo(temp);
		w.close();
		
		r.close();			
	}
	
	@Test
	public void testDistributionMerge() throws IOException {
		File orig = new File(getClass().getResource("/testdata/sampleOneAndTwo.coverage.json").getFile());
		FileReader r = new FileReader(orig);
		IntegerDistributionSet dset = IntegerDistributionSet.read(r);
		IntegerDistributionSet doubleSet = IntegerDistributionSet.merge(dset,dset);
		IntegerDistribution dsetAllDS = dset.getAllSamplesDistributionSet().getAllSequencesDistributionSet().getAllCategoriesDistribution();
		IntegerDistribution doubleSetAllDS = doubleSet.getAllSamplesDistributionSet().getAllSequencesDistributionSet().getAllCategoriesDistribution();
		assertEquals(dsetAllDS.mean(),doubleSetAllDS.mean(),0.0001);
		assertEquals(dsetAllDS.count() * 2 ,doubleSetAllDS.count());
	}
	

	

	
}
