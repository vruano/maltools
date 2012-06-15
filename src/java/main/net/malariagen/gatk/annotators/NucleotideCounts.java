package net.malariagen.gatk.annotators;

import java.util.List;

import net.malariagen.gatk.math.Beta;
import net.malariagen.gatk.math.Phred;
import net.malariagen.utils.Nucleotide;
import net.malariagen.utils.NucleotideIUPAC;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.Allele;

/**
 * Special nucleotide optimized counter.
 * @author valentin
 *
 */
public class NucleotideCounts {
	
	int a;
	int t;
	int g;
	int c;
	int total;
	double errorRateSum;
	double errorRateAvg = Double.NaN;
	int aQsum = 0;
	int tQsum = 0;
	int gQsum = 0;
	int cQsum = 0;
	
	double aQ = Double.NaN;
	double tQ = Double.NaN;
	double gQ = Double.NaN;
	double cQ = Double.NaN;
	
	
	public NucleotideCounts() {
		
	}
	
	public int getCount(Nucleotide x) {
		switch (x) {
		case A: return a;
		case T: return t;
		case C: return c;
		case G: return g;
		default:
			throw new IllegalArgumentException();
		}
	}
	
	public int getScoreSum(Nucleotide x) {
		switch (x) {
		case A: return aQsum;
		case T: return tQsum;
		case C: return cQsum;
		case G: return gQsum;
		default:
			throw new IllegalArgumentException();
		}
	}
	
	public List<Nucleotide> present;

	
	
	public double getQuality(Nucleotide x) {
		switch (x) {
		case A: return aQ;
		case T: return tQ;
		case C: return cQ;
		case G: return gQ;
		default:
			throw new IllegalArgumentException();
		}
	}	
	
	public int getCount(Allele a) {
		return getCount(a.getBases()[0]);
	}
	
	public int getScoreSum(Allele a) {
		return getScoreSum(a.getBases()[0]);
	}
	
	public int getScoreSum(byte b) {
		return getScoreSum(Nucleotide.fromByte(b));
	}
	
	public int getCount(byte b) {
		return getCount(Nucleotide.fromByte(b));
	}
	
	public double getQuality(Allele a) {
		return getQuality(a.getBases()[0]);
	}
	
	public double getQuality(byte b) {
		return getQuality(Nucleotide.fromByte(b));
	}
	


	public void add(AlignmentContext ac) {
		for (PileupElement e : ac.getBasePileup()) {
			NucleotideIUPAC code = NucleotideIUPAC.fromByte(e.getBase());
			byte qual = baseQual(e);
			switch (code) {
			case A: a++; aQsum += qual;  break;
			case U: case T: t++; tQsum += qual; break;
			case G: g++; gQsum += qual; break;
			case C: c++; cQsum += qual;
			default :
			}
			total++;
			errorRateSum += Phred.prob(qual);
		}
	}
	
	public void add(AlignmentContext ...acs) {
		for (AlignmentContext ac : acs)
			add(ac);		
	}
	
	public void add(Iterable<? extends AlignmentContext> acs) {
		for (AlignmentContext ac : acs)
			add(ac);
	}
	
	public void updateQualities() {
		errorRateAvg = errorRateSum / total;
		aQ = alleleQuality(errorRateAvg,total,a);
		tQ = alleleQuality(errorRateAvg,total,t);
		cQ = alleleQuality(errorRateAvg,total,c);
		gQ = alleleQuality(errorRateAvg,total,g);
	}
	
	public void clear() {
		total = 0;
		a = 0;
		t = 0;
		c = 0;
		g = 0;
		errorRateSum = 0;
		errorRateAvg = Double.NaN;
		aQ = tQ = cQ = gQ = Double.NaN;
	}
	
	private double alleleQuality (double er, int total, int count) {	
	    if (count == 0) 
	       return Beta.phred(er, 1.0000001,total + 0.000001);
	    else if (total == count)
	       return Beta.phred(1 - er, 1.000001,count + 0.000001);
	    else
	       return Beta.phred(er,count + 1.000001,total - count + 0.000001);
	    
	}
	
	
	private byte baseQual(PileupElement e) {
		int mq = e.getMappingQual();
		byte bq = e.getQual();
		byte result = mq < bq ? (byte) mq : bq;
		// 3 equal to 0.5 chances. 
		return result < 3 ? 3 : result;
	}
	
	
}