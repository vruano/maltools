package net.malariagen.gatk.genotyper.models.test;

import static org.junit.Assert.*;

import java.util.HashMap;
import java.util.Map;

import net.malariagen.gatk.genotyper.SnpGenotypingContext;
import net.malariagen.gatk.genotyper.models.BernoulliGenotypingModel;
import net.malariagen.gatk.genotyper.models.GenotypingModel;
import net.malariagen.gatk.genotyper.models.GenotypingModelUtils;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.junit.Test;

public class GenotypingModelUtilsTest {
	
	public static final Map<String, GenotypingModel> GOOD_MODELS;
	
	static {
		GOOD_MODELS = new HashMap<String,GenotypingModel>();
		GOOD_MODELS.put("Bernoulli(0.3)",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.3));
		GOOD_MODELS.put("Bernoulli(rpp=0.45)",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.45));
		GOOD_MODELS.put("Bernoulli",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.50000001));
		GOOD_MODELS.put("Bernoulli()",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.50000001));	
		GOOD_MODELS.put("Bernoulli/0.3/",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.3));
		GOOD_MODELS.put("Bernoulli / rpp = 0.45 /",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.45));
		GOOD_MODELS.put("Bernoulli//",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.50000001));
		GOOD_MODELS.put("Bernoulli /  /",new BernoulliGenotypingModel(new SnpGenotypingContext(Allele.create("A",true),Allele.create("T",false)),0.50000001));	
	}
	
	@Test
	public void testGoodModels() {
		for (Map.Entry<String, GenotypingModel> e : GOOD_MODELS.entrySet()) {
			String spec = e.getKey();
			GenotypingModel model = GenotypingModelUtils.getModelInstance(spec);
			assertNotNull("instantiating '" + spec + "'",model);
			assertEquals("comparing '" + spec + "'",e.getValue(),model);
		}
	}

}
