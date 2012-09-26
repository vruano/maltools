package net.malariagen.gatk.uniqueness;

import org.broad.tribble.AbstractFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;

public class UQNCodec extends AbstractFeatureCodec<UQNFeature> {

	@Override
	public UQNFeature decode(String line) {
		if (line.startsWith("#") || line.trim().length() == 0) 
			return null;
		String[] fields = line.split("\\t");
		if (fields.length < 3) 
			throw new UQNInvalidFeatureLineException(line, "not enough fields in line");
		UQNFeature feature = new UQNFeature();
		int nextIndex = 0;
		feature.setSequence(fields[nextIndex++]);
		feature.setPosition(Integer.parseInt(fields[nextIndex++]));
		feature.setScore(Integer.parseInt(fields[nextIndex++]));
		return feature;
	}

	@Override
	public Feature decodeLoc(String line) {
		return decode(line);
	}

	@Override
	public Class<UQNFeature> getFeatureType() {
		return UQNFeature.class;
	}

	@Override
	public Object readHeader(LineReader arg0) {
		return null;
	}

	
	@Override
	public boolean canDecode(String path) {
		if (path.endsWith("uq"))
			return true;
		else 
			return false;
	}
}
