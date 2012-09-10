package net.malariagen.gatk.math;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.List;

import net.malariagen.gatk.coverage.CoverageDistributionException;
import net.malariagen.gatk.math.IntegerDistributionSet;

import org.broadinstitute.sting.commandline.Gatherer;
import org.broadinstitute.sting.utils.exceptions.UserException;

public class IntegerDistributionSetGatherer extends Gatherer {

	@Override
	public void gather(List<File> infiles, File outfile) {
		IntegerDistributionSet[] insets = new IntegerDistributionSet[infiles
				.size()];
		int nextIndex = 0;
		for (File f : infiles)
			try {
				Reader r = new FileReader(f);
				insets[nextIndex++] = IntegerDistributionSet.read(r);
				r.close();
			} catch (FileNotFoundException e) {
				throw new UserException(
						"Input coverage distribution set file '" + f
								+ "' could not be open for reading", e);

			} catch (IOException e) {
				throw new CoverageDistributionException(
						"Error when reading coverage distribution set file", e);
			}
		IntegerDistributionSet merged = IntegerDistributionSet.merge(insets);
		try {
			Writer writer = new FileWriter(outfile);
			merged.write(writer);
			writer.close();
		} catch (IOException e) {
			throw new CoverageDistributionException(
					"Error attempting to write a distribution into a file '"
							+ outfile + "'", e);
		}
	}

	public void gatherTwoAtATime(List<File> inputs, File outfile) {

		try {

			IntegerDistributionSet result = IntegerDistributionSet
					.read(new FileReader(inputs.get(0)));

			for (int i = 1; i < inputs.size(); i++) {
				IntegerDistributionSet other = IntegerDistributionSet
						.read(new FileReader(inputs.get(i)));
				result = IntegerDistributionSet.merge(result, other);
			}
			try {
				Writer writer = new FileWriter(outfile);
				result.write(writer);
				writer.close();
			} catch (IOException e) {
				throw new CoverageDistributionException(
						"Error attempting to write a distribution into a file '"
								+ outfile + "'", e);
			}

		} catch (FileNotFoundException e) {
			throw new UserException("Input coverage distribution set file could not be open for reading", e);

		} catch (IOException e) {
			throw new CoverageDistributionException(
					"Error when reading coverage distribution set file", e);
		}
	}

}
