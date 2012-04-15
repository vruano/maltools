package net.malariagen.gatk.csl;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;

public class CSLCodec implements FeatureCodec<CSLFeature> {

	@Override
	public CSLFeature decode(String line) {
		if (line.startsWith("#") || line.trim().length() == 0) 
			return null;
		String[] fields = line.split("\\t");
		if (fields.length < 4) 
			throw new CSLInvalidFeatureLineException(line, "not enough fields in line");
		if (fields.length > 5) 
			throw new CSLInvalidFeatureLineException(line, "too many fields in line");
		if (fields.length == 5 && !fields[0].equals("SNPS"))
			return null;
		CSLFeature feature = new CSLFeature();
		int nextIndex = fields.length == 5 ? 1 : 0;
		feature.setSequence(fields[nextIndex++]);
		feature.setPosition(Integer.parseInt(fields[nextIndex++]));
		feature.setReference((byte) fields[nextIndex++].charAt(0));
		feature.setAlternative((byte) fields[nextIndex++].charAt(0));
		return feature;
	}

	@Override
	public Feature decodeLoc(String line) {
		return decode(line);
	}

	@Override
	public Class<CSLFeature> getFeatureType() {
		return CSLFeature.class;
	}

	@Override
	public Object readHeader(LineReader arg0) {
		return null;
	}

}
