package net.malariagen.gatk.math;

import static org.junit.Assert.*;

import net.malariagen.gatk.genotyper.models.BernoulliGenotypingModel;

import org.apache.commons.math.MathException;
import org.junit.Test;

public class BetaTest {

	@Test
	public void testLogDoubleDouble() {
		assertEquals(0,Beta.log(1,1),0.0000000001);
		assertEquals(-7.927324,Beta.log(6,6),0.000001);
		assertEquals(9.30565,Beta.log(0.001,0.0001), 0.00001);
		assertEquals(-57857.05,Beta.log(7000,10000121),0.01);
		assertEquals(11.5128,Beta.log(121221,0.00001),0.0001);
	}

	@Test
	public void testLog10DoubleDouble() {
		assertEquals(0,Beta.log10(1,1),0.0000000001);
		assertEquals(-7.927324/Math.log(10),Beta.log10(6,6),0.000001);
		assertEquals(9.30565/Math.log(10),Beta.log10(0.001,0.0001), 0.00001);
		assertEquals(-57857.05/Math.log(10),Beta.log10(7000,10000121),0.01);
		assertEquals(11.5128/Math.log(10),Beta.log10(121221,0.00001),0.0001);
	}

	@Test
	public void testLogDoubleDoubleDouble() throws MathException {
		assertEquals(0,Beta.log(0,1,1),0.0000000001);
		assertEquals(1,Beta.log(1,1,1),0.0000000001);
		assertEquals(-.3566749,Beta.log(0.7,1,1),0.0000001);
        assertEquals(-.09431068,Beta.log(0.7, 1,2),0.0000001);
        assertEquals(-.000100005,Beta.log(0.99,1,2),0.000001);
        assertEquals(-1e-08,Beta.log(0.9999, 1,2),0.000001);
        assertEquals(-10.81978,Beta.log(0.00001, 1,2),0.00001);
        assertEquals(-7.133499,Beta.log(0.7, 20,1),0.000001);
        assertEquals(-.2010067,Beta.log(0.99,20,1),0.000001);
        assertEquals(-0.00200001,Beta.log(0.9999, 20,1),0.000001);
        assertEquals(-230.2585,Beta.log(0.00001, 20,1),0.0001);
	}

	@Test
	public void testLog10DoubleDoubleDouble() {
		assertEquals(0,Beta.log10(0,1,1),0.0000000001);
		assertEquals(1,Beta.log10(1,1,1),0.0000000001);
		assertEquals(-10.97384,Beta.log10(0.0036014098530065546,6,10),0.00001);
		assertEquals(-4.612566e-12,Beta.log10(1-0.00360141,10,6),1.e-18);
		assertEquals(-.3566749/Math.log(10),Beta.log10(0.7,1,1),0.0000001);
        assertEquals(-.09431068/Math.log(10),Beta.log10(0.7, 1,2),0.0000001);
        assertEquals(-.000100005/Math.log(10),Beta.log10(0.99,1,2),0.000001);
        assertEquals(-1e-08/Math.log(10),Beta.log10(0.9999, 1,2),0.000001);
        assertEquals(-10.81978/Math.log(10),Beta.log10(0.00001, 1,2),0.00001);
        assertEquals(-7.133499/Math.log(10),Beta.log10(0.7, 20,1),0.000001);
        assertEquals(-.2010067/Math.log(10),Beta.log10(0.99,20,1),0.000001);
        assertEquals(-0.00200001/Math.log(10),Beta.log10(0.9999, 20,1),0.000001);
        assertEquals(-230.2585/Math.log(10),Beta.log10(0.00001, 20,1),0.0001);

	}

}
