package net.malariagen.vcf;

import java.util.Collections;
import java.util.Map;

import org.broadinstitute.sting.utils.codecs.vcf.VCFSimpleHeaderLine;

public class VCFReadGroupHeaderLine  extends VCFSimpleHeaderLine {

	public VCFReadGroupHeaderLine(String name,
			Map<String, String> genericFields) {
		super("readGroup", name, genericFields);
	}

	@SuppressWarnings("unchecked")
	public VCFReadGroupHeaderLine(String name) {
		super("readGroup", name,  Collections.EMPTY_MAP);
	}


}
