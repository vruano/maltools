package net.malariagen.gatk.genotyper.models.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Formatter;
import java.util.List;

import net.malariagen.gatk.genotyper.SnpGenotypingContext;
import net.malariagen.gatk.genotyper.models.DiscreteMixtureGenotypingModel;
import net.malariagen.gatk.test.WalkerTest;

import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.junit.Test;

public class DiscreteMixturePolyploidModelTest extends WalkerTest {

    public static final String PRIOR_1PC_RESOURCE = "mockups/discrete-mixture/prior.1pc";
    public static final String PRIOR_1PC_SIMPLE_RESOURCE = "mockups/discrete-mixture/prior.1pc.simple";    
    public static final String PRIOR_10PC_RESOURCE = "mockups/discrete-mixture/prior.10pc";
    public static final String PRIOR_5PC_RESOURCE = "mockups/discrete-mixture/prior.5pc";


	@Test
	public void testPriorLoad() {
		DiscreteMixtureGenotypingModel dmgm1pc = modelFromPrior(PRIOR_1PC_RESOURCE);
		DiscreteMixtureGenotypingModel dmgm10pc = modelFromPrior(PRIOR_10PC_RESOURCE);
		DiscreteMixtureGenotypingModel dmgm5pc = modelFromPrior(PRIOR_5PC_RESOURCE);
	}


	private DiscreteMixtureGenotypingModel modelFromPrior(String res) {
		File prior1pc = new File(getClass().getClassLoader().getResource(res).getFile());
		DiscreteMixtureGenotypingModel dmgm1pc = new DiscreteMixtureGenotypingModel(prior1pc,0);
		dmgm1pc.setGenotypingContext(new SnpGenotypingContext(Allele.create((byte) 'A',true),Allele.create((byte)'T',false)));
		return dmgm1pc;
	}
	
	
	@Test
	public void testModelInstanciation() {
		DiscreteMixtureGenotypingModel dmgm = new DiscreteMixtureGenotypingModel();
	}
	
	@Test
	public void test7G8xGb4() throws IOException {
		File prior1pc = new File(getClass().getClassLoader().getResource(PRIOR_1PC_RESOURCE).getFile());
		File baseDir = new File("/data/haldane/malariagen/pfalciparum/crosses/7g8xGb4/som/");
		File _7g8 = new File(baseDir, "PG0084-C.bam");
		File _Gb4 = new File(baseDir, "PG0083-C.bam");

		String[] progenyNames = { "PG0100-C", "PG0101-C", "PG0092-C",
				"PG0112-C", "PG0113-C", "PG0110-C", "PG0111-C", "PG0085-C",
				"PG0095-C", "PG0103-C", "PG0094-C", "PG0087-C", "PG0104-C",
				"PG0107-C", "PG0096-C", "PG0098-C", "PG0106-C", "PG0109-C",
				"PG0088-C", "PG0086-C", "PG0099-C", "PG0093-C", "PG0091-C",
				"PG0105-C", "PG0090-C", "PG0097-C", "PG0108-C", "PG0102-C" };
		File[] bamFiles = new File[progenyNames.length + 2];
		bamFiles[0] = _7g8;
		bamFiles[1] = _Gb4;
		for (int i = 2; i < bamFiles.length; i++)
			bamFiles[i] = new File(baseDir, progenyNames[i - 2] + ".bam");

		StringBuffer sb = new StringBuffer(1000);

		for (File f : bamFiles)
			sb.append("-I ").append(f.toString()).append(' ');

		File reference = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference/3D7_pm.fa").getFile());
		File uniqness = new File(getClass().getResource(
				"/mockups/bernoulli-test/reference/3D7_pm.uq").getFile());
		File output = File.createTempFile("BernoulliTest", ".vcf");
		Formatter formatter = new Formatter();
		Formatter cmdSpec = formatter
				.format(" -R %s -T MetaGenotyper -o %s -baseq_do %s.bdo -mbq 20 -mmq 20 -out_mode EMIT_ALL_SITES -gem EMIT_ALL -smodel DiscreteMixture[prior=%s] -A AverageBaseQuality -A ReadDepthAndAllelicFractionBySample -A DepthPerAlleleByVariant -A UniquenessScore -A NumberSamplesWithData -B:uniqueness,UQN %s %s",
						reference, output, output, prior1pc.toString(), uniqness, sb.toString());
		// No md5s for now.
		List<String> md5s = Collections.emptyList();
		WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
				cmdSpec.toString(), md5s);
		executeTest("testCoverageCounting", spec);
	}


}
