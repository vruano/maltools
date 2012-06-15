package net.malariagen.gatk.genotyper;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import net.malariagen.gatk.annotators.NucleotideCounts;
import net.malariagen.utils.NucleotideIUPAC;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.Allele;

public class SnpGenotypingContext implements GenotypingContext {


	private Allele referenceAllele;

	private List<Allele> allAlleles;
	
	private NucleotideCounts counts;
	

	public SnpGenotypingContext(ReferenceContext rc,
			AlignmentContext... acs) {
		counts = new NucleotideCounts();
		counts.add(acs);
		
		referenceAllele = Allele.create(rc.getBase(),true);
		allAlleles = new ArrayList<Allele>(4);
		allAlleles.add(referenceAllele);
		
		byte[] protoBytes = new byte[10]; // rarely is going to be more than
		int[] altCounts = new int[10];

		if (count > 2 && biallelic) {
			int maxAltCount = altCounts[0];
			int maxAltIdx = 0;
			for (int i = 1; i < count - 1; i++)
				if (altCounts[i] > maxAltCount)
					maxAltCount = altCounts[maxAltIdx = i];
			protoBytes[1] = protoBytes[maxAltIdx + 1];
			count = 2;
			protoBytes = Arrays.copyOf(protoBytes, 2);
		}
		Arrays.sort(protoBytes, 1, count);
		bytes = protoBytes;
		
	}

	

	private byte[] bytes;

	private int hashCode;


	public SnpGenotypingContext(boolean biallelic, ReferenceContext rc,
			AlignmentContext[] acs, byte possibleAlternatives) {
		byte[] protoBytes = new byte[10]; // rarely is going to be more than
											// 10.
		int[] altCounts = new int[10];
		int count = 0;
		protoBytes[count++] = rc.getBase();
		for (AlignmentContext ac : acs) {
			for (PileupElement e : ac.getBasePileup()) {
				byte b = e.getBase();
				if (!NucleotideIUPAC.areCompatible(b, possibleAlternatives))
					continue;
				boolean found = false;
				for (int i = 0; i < count; i++)
					if (protoBytes[i] == b) {
						if (i > 0)
							altCounts[i - 1] += e.getQual();
						found = true;
						break;
					}
				if (!found) {
					if (count == protoBytes.length) // Just in case!!
						protoBytes = Arrays.copyOf(protoBytes,
								protoBytes.length << 1);
					altCounts[count - 1] = e.getQual();
					protoBytes[count++] = b;
				}
			}
		}
		if (count > 2 && biallelic) {
			int maxAltCount = altCounts[0];
			int maxAltIdx = 0;
			for (int i = 1; i < count - 1; i++)
				if (altCounts[i] > maxAltCount)
					maxAltCount = altCounts[maxAltIdx = i];
			protoBytes[1] = protoBytes[maxAltIdx + 1];
			count = 2;
			protoBytes = Arrays.copyOf(protoBytes, 2);
		}
		Arrays.sort(protoBytes, 1, count);
		bytes = protoBytes;
	}

	public int hashCode() {
		return hashCode;
	}

	public boolean equals(Object o) {
		if (o instanceof byte[])
			return Arrays.equals(bytes, (byte[]) o);
		else
			return false;
	}
	
	
	public boolean hasAllele(Allele a) {
		return hasAllele(a.getBases()[0]);
	}

	public boolean hasAllele(byte b) {
		if (b == 'n')
			return false;
		if (b == 'N')
			return false;
		if (b == '.')
			return false;
		if (b == '-')
			return false;
		for (int i = 0; i < allAlleles.size(); i++) {
			Allele a2 = allAlleles.get(i);
			byte b2 = a2.getBases()[0];
			if (b == b2)
				return true;
		}
		return false;
	}




	@Override
	public Allele getReferenceAllele() {
		return referenceAllele;
	}

	@Override
	public int getAlleleCount() {
		return allAlleles.size();
	}

	@Override
	public Allele getAllele(int index) {
		return allAlleles.get(index);
	}

	@Override
	public Allele getAllele(byte b) {
		byte[] ba = bytes;
		if (b == ba[0])
			return referenceAllele;
		else {
			for (int i = 0; i < ba.length; i++)
				if (b == ba[i])
					return allAlleles.get(i);
			return null;
		}
	}

	@Override
	public int getAlleleIndex(byte b) {
		byte[] ba = bytes;
		for (int i = 0; i < ba.length; i++)
			if (ba[i] == b)
				return i;
		return -1;
	}

	@Override
	public List<Allele> getAlleleList() {
		return allAlleles;
	}

	private static class Key {

	}

}
