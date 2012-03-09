package net.malariagen.gatk.math;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;

public class CombineIntegerDistributionSets extends CommandLineProgram {

    @Input(fullName = "inputs", shortName = "I", doc = "Input integer distribution JSON files to merge", required = true)
    public List<File> inputs = new ArrayList<File>();

    @Output(fullName = "output", shortName = "o", doc = "Output integer distribution JSON file", required = true)
    public File output;
        
	@Override
	protected int execute() throws Exception {
		IntegerDistributionSetGatherer gatherer = new IntegerDistributionSetGatherer();
		gatherer.gatherTwoAtATime(inputs, output);
		return 0;
	}

	public static int executeWitoutExiting(String[] args) throws Exception {
			CombineIntegerDistributionSets instance = new CombineIntegerDistributionSets();
			start(instance,args);
		return CommandLineProgram.result;
	}	
	
	public static void main(String[] args) {
		try {
			CombineIntegerDistributionSets instance = new CombineIntegerDistributionSets();
			start(instance,args);
		} catch (Exception e) {
			exitSystemWithError(e);
		}
		System.exit(CommandLineProgram.result);
	}

}
