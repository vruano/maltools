package net.malariagen.gatk.genotyper;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import net.malariagen.gatk.annotators.CoverageAnnotation;
import net.malariagen.gatk.math.IntegerDistributionSet;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Custom diploid genotype caller based on the {@link UnifiedGenotyper}.
 * 
 * For now it only adds the possibility of annotate coverage depth where the distribution is passed as an argument (-CvgD or -coverageDistribution).
 * 
 * @author valentin
 *
 */
public class CoverageAnnotationEngine {

	
	@Input(fullName = "coverageDistribution", shortName = "CvgD", doc = "File where to find the precalcuated coverage distribution", required = false)
	public File coverageDistributionFile;
	
	
	public CoverageAnnotationEngine(List<String> annotationToUse, File coverageDistributionFile) {
		this.annotationsToUse = annotationToUse;
		this.coverageDistributionFile = coverageDistributionFile;
	}

	private List<String> annotationsToUse;
	private IntegerDistributionSet coverageDistribution;

	public void initializeCoverageDistribution() {
		if (coverageDistributionFile == null) {
			for (String n : CoverageAnnotation.ANNOTATION_NAMES) {
				if (annotationsToUse.contains(n)) {
					Logger.getLogger(Walker.class).warn("no coverage distribution file provided despite annotation requested. Coverage won't be annotated");
					annotationsToUse
							.removeAll(CoverageAnnotation.ANNOTATION_NAMES);
					break;
				}
			}
		} else {
			if (!coverageDistributionFile.exists())
				new UserException("the coverage distribution file provided '"
						+ coverageDistributionFile.toString()
						+ "' does not seem to exists");
			if (!coverageDistributionFile.isFile())
				new UserException("the coverage distribution file provided '"
						+ coverageDistributionFile.toString()
						+ "' does not seem to be a regular file");
			if (!coverageDistributionFile.canRead())
				new UserException("the coverage distribution file provided '"
						+ coverageDistributionFile.toString()
						+ "' cannot be read");
			try {
				coverageDistribution = IntegerDistributionSet
						.read(new FileReader(coverageDistributionFile));
			} catch (IOException e) {
				throw new UserException(
						"the coverage distribution file provided '"
								+ coverageDistributionFile.toString()
								+ "' could not be read to completion ", e);
			}

			// Add the corresponding annotation if not already provided by the
			// user.
			boolean annotationFound = false;
			for (String n : CoverageAnnotation.ANNOTATION_NAMES) {
				if (annotationsToUse.contains(n)) {
					annotationFound = true;
					break;
				}
			}
			if (!annotationFound) 
				annotationsToUse.add(CoverageAnnotation.class.getSimpleName());
		}
	}
	
	public void setupCoverageAnnotation() {
		if (coverageDistribution != null)
			CoverageAnnotation.coverageDistributionSet.set(coverageDistribution);
	}

}
