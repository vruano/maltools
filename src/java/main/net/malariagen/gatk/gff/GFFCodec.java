package net.malariagen.gatk.gff;

import org.broad.tribble.AbstractFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineReader;

public class GFFCodec extends AbstractFeatureCodec<GFFFeature> {

	@Override
	public Feature decodeLoc(String line) {
		return decode(line);
	}

	@Override
	public GFFFeature decode(String line) {
		if (line.startsWith("#") || line.trim().length() == 0) 
			return null;
		String[] fields = line.split("\\t");
		if (fields.length < 8)
			throw new GFFInvalidFeatureLineException(line, "too few fields");
		int nextIndex = 0;
		String sequenceName = fields[nextIndex++];
		String source = fields[nextIndex++];
		String typeString = fields[nextIndex++];
		int start;
		int end;
		try {
			start = Integer.parseInt(fields[nextIndex++]);
			end = Integer.parseInt(fields[nextIndex++]);
		} catch (NumberFormatException e) {
			throw new GFFInvalidFeatureLineException(line,
					"invalid numeric value");
		}
		if (start > end)
			throw new GFFInvalidFeatureLineException(line,
					"feature start cannot be greater than the end");
		String scoreString = fields[nextIndex++];
		GFFStrand strand = GFFStrand.fromCharSequence(fields[nextIndex++]);
		GFFFrame frame = GFFFrame.fromCharSequence(fields[nextIndex++]);
		String attributeString = fields.length > nextIndex ? fields[nextIndex++]
				: "";

		GFFFeature feature = new GFFFeature();
		feature.setSequenceName(sequenceName);
		feature.setStart(start);
		feature.setEnd(end);
		feature.setTypeString(typeString);
		feature.setSource(source);
		feature.setFrame(frame);
		try {
			feature.setScore(scoreString);
		} catch (NumberFormatException e) {
			throw new GFFInvalidFeatureLineException(line,
					"invalid score value");
		}
		feature.setStrand(strand);
		feature.setAttributeString(attributeString);

		return feature;
	}

	@Override
	public Class<GFFFeature> getFeatureType() {
		return GFFFeature.class;
	}

	@Override
	public Object readHeader(LineReader reader) {
		return null;
	}
	
	@Override
	public boolean canDecode(String path) {
		if (path.endsWith("gff"))
			return true;
		else 
			return false;
	}

}
