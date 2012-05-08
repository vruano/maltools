package Maltools::Tools::Gatk::SgGenotyper;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) }
  
our $INPUTS = {
    ref => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
    ref_uq => { type => 'file' , mandatory => 1 },
    ref_gff => { type => 'file', mandatory => 1 },
    samples => { type => { name => 'bam', indexed => 1}, multiple => 1, mandatory => 1 },
    regions => { type => 'string', multiple => 1, mandatory => 0 },
    interval_list => { type => 'file' },
    threads => { type => 'num' ,default => 100 },
};

our $OUTPUTS = { 
    out => { type => 'vcf', mandatory => 1 },
    filtered_out => { type => 'vcf' },
    out_eval => { type => 'file' },
    filtered_out_eval => { type => 'file' },
    indel_out => { type => 'vcf' },
    cvg_distribution => { type => 'file', mandatory => 1 },
};


sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub job {
  my ($self,%args) = @_;
  $args{cpu_ratio} = $self->calculate_cpu_ratio(%args) unless exists $args{cpu_ratio};
  $args{memory} = $self->calculate_memory(%args) unless exists $args{memory};
  $args{cpu_time} = $self->calculate_cpu_time(%args) unless exists $args{cpu_time};
  return $self->SUPER::job(%args);
}

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 0.1;
}

sub calculate_memory {
  my ($self,%args) = @_;
  return 2000;
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  return 3 * 60 * 60;
}

sub interpreter {
   my ($self) = @_;
   return 'queue';
}

sub command_template {
  my ($self,%args) = @_;
  TTS('{$J->interpreter()} --memory {$J->memory} -S {$S} ' 
     . '-R {$J->input("ref")} '
     . '-A DepthOfCoverage '
     . '-A AbsoluteCounts '
     . '-A UniquenessScore '
     . '-A CodingAnnotation '
     . '-A CoverageAnnotation '
     . '-cvgO {$J->output("cvg_distribution")} ' 
     . '{$il = $J->input("interval_list"); $il ? "-L $il " : " " }'
     . '-filter LowCoverage '
     . "-filterExpression \"'CMF < 0.5'\" "
     . '-filter HighCoverage '
     . "-filterExpression \"'CMF > 2'\" "
     . '-filter HARD_TO_VALIDATE '
     . "-filterExpression \"'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)'\" "
     . '-filter LowMrAF '
     . "-filterExpression \"'MrAF < 0.01'\" "
     . '-filter LowMrAmaxSD '
     . "-filterExpression \"'MrAmaxSD < 10'\" "
     . ' { join(" ",map {"-I $_"} @{$J->input("samples")} ); } '
     . ' { join(" ",map {"-L $_"} @{$J->input("regions")} ); } '
     . '-o {$J->output("out")} '
     . '-jobQueue long '
     . '-bsub -run' );
}


1;


__DATA__


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An example building on the intro ExampleCountReads.scala.
 * Runs an INCOMPLETE version of the UnifiedGenotyper with VariantEval and optional VariantFiltration.
 */
class ExampleUnifiedGenotyper extends QScript \{
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFiles: List[File] = Nil

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Output(doc="output file",shortName="o", required=false)
  var output: File = _

  @Output(fullName="coverageFile", doc="coverage output file",shortName="cvgO", required=true)
  var coverageFile: File = _ 

  @Argument(doc="A optional list of filter names.", shortName="filter", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(doc="An optional list of filter expressions.", shortName="filterExpression", required=false)
  var filterExpressions: List[String] = Nil

  @Argument(doc="List of additional annotation.", shortName="A", required=false)
  var annotations: List[String] = Nil 

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK \{
    this.reference_sequence = qscript.referenceFile
    this.intervals = List(qscript.intervals)
    //  is how you set the value for an scala Option.
    // Set the memory limit to 2 gigabytes on each command.
    this.memoryLimit = 2000
  \}


  def script = \{
    // Create the four function that we can run.
    val coverageCounter = new CountCoverage with UnifiedGenotyperArguments
    val genotyper = new DiploidGenotyper with UnifiedGenotyperArguments
    val indelGenotyper = new UnifiedGenotyper with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    coverageCounter.scatterCount = 32
    coverageCounter.input_file = qscript.bamFiles
    coverageCounter.out = qscript.coverageFile
    coverageCounter.countSamples = true
    coverageCounter.rodBind :+= RodBind("features","GFF", new File("{$J->input("ref_gff")}"))

    indelGenotyper.scatterCount = {$J->input("threads")}
    indelGenotyper.input_file = qscript.bamFiles
    indelGenotyper.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    indelGenotyper.out = new File("{$J->output("indel_out")}")
    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    genotyper.scatterCount = {$J->input("threads")}
    genotyper.input_file = qscript.bamFiles
    genotyper.stand_emit_conf = 3.0
    genotyper.coverageDistribution = coverageCounter.out
    genotyper.annotation = qscript.annotations
    genotyper.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP

    genotyper.rodBind :+= RodBind("uniqueness", "UQN", new File("{$J->input("ref_uq")}"))
    genotyper.rodBind :+= RodBind("features", "GFF", new File("{$J->input("ref_gff")}"))
    
    genotyper.out = new File("{$J->output("out")}")

    // qscript.output 

    
    evalUnfiltered.rodBind :+= RodBind("eval", "VCF", genotyper.out)
    evalUnfiltered.out = new File("{$J->output("out_eval")}") 
    // swapExt(genotyper.out, "vcf", "eval")

    variantFilter.rodBind :+= RodBind("variant", "VCF", genotyper.out)
    variantFilter.rodBind :+= RodBind("mask", "VCF", indelGenotyper.out)
    variantFilter.maskName = "AroundIndel"
    variantFilter.clusterWindowSize = 10

    variantFilter.out = new File("{$J->output("filtered_out")}")
    // swapExt(qscript.output, "vcf", "filtered.vcf")
    variantFilter.filterName = filterNames
    variantFilter.filterExpression = filterExpressions

    evalFiltered.rodBind :+= RodBind("eval", "VCF", variantFilter.out)
    evalFiltered.out = new File("{$J->output("filtered_out_eval")}") 
    // swapExt(variantFilter.out, "vcf", "eval")
    //add(indelGenotyper,coverageCounter);
    add(coverageCounter);
    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    //if (filterNames.size > 0)
    add(variantFilter, evalFiltered)
  \}
\}
