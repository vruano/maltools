package net.malariagen.gatk.gff;

import java.util.HashMap;
import java.util.Map;

public enum GFFType {
	ATTENUATOR, C_REGION("C_region"), CAAT_SIGNAL("CAAT_signal"), CDS("CDS"), D_LOOP(
			"D-loop"), D_SEGMENT("D_segment"), ENHANCER, EXON, GAP, GC_SIGNAL(
			"GC_signal"), GENE, I_DNA("iDNA"), INTRON, J_SEGMENT("J_segment"), LTR(
			"LTR"), MAT_PEPTIDE, MISC_BINDING, MISC_DIFFERENCE, MISC_FEATURE, MISC_RNA(
			"misc_RNA"), MISC_SIGNAL, MISC_STRUCTURE, MOBILE_ELEMENT, MODIFIED_BASE, M_RNA(
			"mRNA"), NC_RNA("ncRNA"), N_REGION("N_region"), OLD_SEQUENCE, OPERON, ORIT(
			"oriT"), POLYA_SIGNAL("polyA_signal"), POLYA_SITE("polyA_site"), PRECURSOR_RNA(
			"precursor_RNA"), PRIM_TRANSCRIPT, PRIMER_BIND, PROMOTER, PROTEIN_BIND, RBS(
			"RBS"), REPEAT_REGION, REP_ORIGIN, R_RNA("rRNA"), S_REGION(
			"S_region"), SIG_PEPTIDE, SOURCE, STEM_LOOP, STS("STS"), TATA_SIGNAL(
			"TATA_signal"), TERMINATOR, TM_RNA("tmRNA"), TRANSIT_PEPTIDE, T_RNA(
			"tRNA"), UNSURE, V_REGION("V_region"), V_SEGMENT("V_segment"), VARIATION, _3_UTR(
			"3'UTR"), _5_UTR("5'UTR"), _10_SIGNAL("-10_signal"), _35_SIGNAL(
			"-35_signal"), UNKNOWN("<unknown>");

	// indicates the best way to represent this is a string (with non id chars
	// like - or mixed upper-lower case).
	private String preferredString;

	// Constructs a type with the default preferred string representation
	GFFType() {
		this(null);
	}

	// Constructs a feature type with an arbitrary preferred representation.
	// Null indicates that the lower case of the default string representation
	// is appropriate for it.
	GFFType(String preferredString) {
		if (preferredString == null)
			preferredString = super.toString().toLowerCase();
		this.preferredString = preferredString;
	}

	// Map from pf. string to type
	private static Map<String, GFFType> stringToValue = new HashMap<String, GFFType>();

	// Initializes the map from pf. string to type.
	// It cannot be done in the constructor due to language restrictions.
	static {
		for (GFFType ft : GFFType.values()) {
			stringToValue.put(ft.toString().toLowerCase(), ft);
			stringToValue.put(ft.preferredString, ft);
		}
	};

	/**
	 * Parses a character sequence into the corresponding feature type.
	 * 
	 * @param cs
	 *            the character sequence to parse.
	 * @throw NullPointerException if {@code cs} is {@code null}.
	 * @return never {@code null} but perhaps {@link #UNKNOWN} if it cannot find
	 *         a suitable feature.
	 */
	public static GFFType fromCharSequence(CharSequence cs) {
		GFFType ft = stringToValue.get(cs.toString().toLowerCase());
		return ft == null ? UNKNOWN : ft;
	}

	// Returns a preferable human friendly representation of this type
	@Override
	public String toString() {
		return preferredString;
	}

	// A false return does not mean that cannot code for a protein but that it
	// is unknown.
	public boolean isProteinCoding() {
		switch (this) {
		case CDS:
			return true;

		default:
			return false;
		}
	}

}
