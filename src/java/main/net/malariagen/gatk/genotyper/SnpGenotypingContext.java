package net.malariagen.gatk.genotyper;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.Allele;

public class SnpGenotypingContext implements GenotypingContext {

	private final Key key;

	private Allele referenceAllele;

	private List<Allele> allAlleles;

	public SnpGenotypingContext(ReferenceContext rc, AlignmentContext... ac) {
		this(false, rc, ac);
	}

	public SnpGenotypingContext(boolean biallelic, ReferenceContext rc,
			AlignmentContext... ac) {
		this(new Key(biallelic, rc, ac));
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

	public SnpGenotypingContext(Allele ref, Allele... alternative) {
		this(new Key(false, ref, alternative));
	}

	protected SnpGenotypingContext(Key key) {
		if (key.bytes.length < 1)
			throw new IllegalArgumentException("missing reference allele byte");
		referenceAllele = Allele.create(key.bytes[0], true);
		allAlleles = new ArrayList<Allele>(key.bytes.length);
		allAlleles.add(referenceAllele);
		byte[] testByte = new byte[1];
		for (int i = 1; i < key.bytes.length; i++) {
			testByte[0] = key.bytes[i];
			if (Allele.acceptableAlleleBases(testByte)
					&& testByte[0] != (byte) 'n' && testByte[0] != (byte) 'N')
				allAlleles.add(Allele.create(testByte[0], false));
		}
		allAlleles = Collections.unmodifiableList(allAlleles);
		this.key = key;
	}

	public SnpGenotypingContext(boolean b, byte allowed,
			ReferenceContext refContext, AlignmentContext... rawContext) {
		this(new Key(b, refContext, rawContext, allowed));
	}

	@Override
	public int hashCode() {
		return key.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof SnpGenotypingContext)
			return key.equals(((SnpGenotypingContext) o).key);
		else
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
		byte[] ba = key.bytes;
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
		byte[] ba = key.bytes;
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

		private byte[] bytes;

		private int hashCode;

		public Key(boolean biallelic, Allele rc, Allele... acs) {
			this(biallelic, rc, acs, (byte) 'N');
		}

		public Key(boolean biallelic, Allele rc, Allele[] acs,
				byte possibleAlternatives) {
			byte[] protoBytes = new byte[10]; // rarely is going to be more than

			int[] altCounts = new int[10];
			int count = 0;
			protoBytes[count++] = rc.getBases()[0];
			for (Allele ac : acs) {
				byte b = ac.getBases()[0];
				if (!compatibleCode(b, possibleAlternatives))
					continue;
				boolean found = false;
				for (int i = 0; i < count; i++)
					if (protoBytes[i] == b) {
						found = true;
						if (i > 0)
							altCounts[i - 1]++;
						break;
					}
				if (!found) {
					if (count == protoBytes.length) // Just in case!!
						protoBytes = Arrays.copyOf(protoBytes,
								protoBytes.length << 1);
					altCounts[count - 1] = 1;
					protoBytes[count++] = b;
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

		private boolean compatibleCode(byte b, byte allowed) {
			switch (b) {
			case 'a':
			case 'A':
				switch (allowed) {
				case 'A':
				case 'N':
				case 'X':
				case 'M':
				case 'R':
				case 'W':
				case 'V':
				case 'H':
				case 'D':
					return true;
				default:
					return false;
				}
			case 'T':
			case 't':
				switch (allowed) {
				case 'T':
				case 'N':
				case 'X':
				case 'W':
				case 'Y':
				case 'K':
				case 'H':
				case 'D':
				case 'B':
					return true;
				default:
					return false;
				}
			case 'c':
			case 'C':
				switch (allowed) {
				case 'C':
				case 'N':
				case 'X':
				case 'M':
				case 'S':
				case 'Y':
				case 'V':
				case 'H':
				case 'B':
					return true;
				default:
					return false;
				}
			case 'G':
			case 'g':
				switch (allowed) {
				case 'G':
				case 'N':
				case 'X':
				case 'R':
				case 'S':
				case 'K':
				case 'V':
				case 'D':
				case 'B':
					return true;
				default:
					return false;
				}
			default:
				return false;
			}
		}

		public Key(boolean biallelic, ReferenceContext rc,
				AlignmentContext... acs) {
			this(biallelic, rc, acs, (byte) 'N');
		}

		public Key(boolean biallelic, ReferenceContext rc,
				AlignmentContext[] acs, byte possibleAlternatives) {
			byte[] protoBytes = new byte[10]; // rarely is going to be more than
												// 10.
			int[] altCounts = new int[10];
			int count = 0;
			protoBytes[count++] = rc.getBase();
			for (AlignmentContext ac : acs) {
				for (PileupElement e : ac.getBasePileup()) {
					byte b = e.getBase();
					if (!compatibleCode(b, possibleAlternatives))
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

	}

}
