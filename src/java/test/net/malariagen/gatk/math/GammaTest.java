package net.malariagen.gatk.math;

import static org.junit.Assert.*;

import org.junit.Test;

public class GammaTest {

	@Test
	public void testLog() {
		assertEquals(Gamma.log(0),Double.POSITIVE_INFINITY,0.000001);
		assertEquals(Gamma.log(1),0,0.0000000001);
		assertEquals(Gamma.log(10),Math.log(362880),0.00001);
	    assertEquals(Gamma.log(235),1046.192,0.001);
	}
	
	public void testLog10() {
		assertEquals(Gamma.log10(0),Double.POSITIVE_INFINITY,0.000001);
		assertEquals(Gamma.log10(1),0,0.0000000001);
		assertEquals(Gamma.log10(10),Math.log10(362880),0.00001);
	    assertEquals(Gamma.log10(235),454.3555,0.0001);
	}

}
