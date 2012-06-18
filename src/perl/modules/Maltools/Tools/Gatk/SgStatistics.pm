package Maltools::Tools::Gatk::SgStatistics;

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
  bq_out => { type => 'file' },
  baq_bq_out => { type => 'file' },
  rbsq_out => { type => 'file' },
  rmq_out => { type => 'file' },
  cvg_out => { type => 'file' },
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
  return 4000;
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  return 8 * 60 * 60;
}

sub interpreter {
   my ($self) = @_;
   return 'queue';
}

sub command_template {
  my ($self,%args) = @_;
  TTS('{$J->interpreter()} --memory {$J->memory} -S {$S} ' 
     . '-R {$J->input("ref")} '
     . '{ $J->output("bq_out") ? "-bqOut " . $J->output("bq_out") . " " : ""}'
     . '{ $J->output("rbsq_out") ? "-rbsqOut " . $J->output("rbsq_out") . " " : ""}'
     . '{ $J->output("rmq_out") ? "-rmqOut " . $J->output("rmq_out") . " " : ""}'
     . '{ $J->output("baq_bq_out") ? "-baqBqOut " . $J->output("baq_bq_out") . " " : ""}'
     . '{ $J->output("cvg_out") ? "-cvgOut " . $J->output("cvg_out") . " " : ""} ' 
     . '{$il = $J->input("interval_list"); $il ? "-L $il " : " " }'
     . ' { join(" ",map {"-I $_"} @{$J->input("samples")} ); } '
     . ' { join(" ",map {"-L $_"} @{$J->input("regions")} ); } '
     . '-jobQueue long '
     . '-bsub -run' );
}


1;


__DATA__

import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode
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

  @Output(doc="output file",fullName="baseQualityOutput",shortName="bqOut", required=false)
  var baseQualityOutput: File = _

  @Output(doc="output file",fullName="baqBaseQualityOutput",shortName="baqBqOut", required=false)
  var baqBaseQualityOutput: File = _
 
  @Output(doc="output file",fullName="readBaseSumQualityOutput",shortName="rbsqOut", required=false)
  var readBaseSumQualityOutput: File = _

  @Output(doc="output file",fullName="readMappingQualityOutput",shortName="rmqOut", required=false)
  var readMappingQualityOutput: File = _
  
  @Output(fullName="coverageFile", doc="coverage output file",shortName="cvgOut", required=false)
  var coverageOutput: File = _ 

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
    val baseQualityCounter = new CountBaseQuality with UnifiedGenotyperArguments
    val baqBaseQualityCounter = new CountBaseQuality with UnifiedGenotyperArguments
    val readQualityCounter = new CountReadQuality with UnifiedGenotyperArguments

    baseQualityCounter.scatterCount = 64
    baseQualityCounter.input_file = qscript.bamFiles
    baseQualityCounter.out = qscript.baseQualityOutput
    baseQualityCounter.countSamples = true
    baseQualityCounter.rodBind :+= RodBind("features","GFF", new File("{$J->input("ref_gff")}"))

    baqBaseQualityCounter.baq = CalculationMode.CALCULATE_AS_NECESSARY
    baqBaseQualityCounter.scatterCount = 64
    baqBaseQualityCounter.input_file = qscript.bamFiles
    baqBaseQualityCounter.out = qscript.baqBaseQualityOutput
    baqBaseQualityCounter.countSamples = true
    baqBaseQualityCounter.rodBind :+= RodBind("features","GFF", new File("{$J->input("ref_gff")}"))

    readQualityCounter.scatterCount = 64
    readQualityCounter.input_file = qscript.bamFiles
    readQualityCounter.baseQualitySumOutput = qscript.readBaseSumQualityOutput
    readQualityCounter.mappingQualityOutput = qscript.readMappingQualityOutput
    readQualityCounter.countSamples = true
    readQualityCounter.rodBind :+= RodBind("features","GFF", new File("{$J->input("ref_gff")}"))

    coverageCounter.scatterCount = 64
    coverageCounter.input_file = qscript.bamFiles
    coverageCounter.out = qscript.coverageOutput
    coverageCounter.countSamples = true
    coverageCounter.rodBind :+= RodBind("features","GFF", new File("{$J->input("ref_gff")}"))

    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.

    if (qscript.baqBaseQualityOutput != null)  add(baqBaseQualityCounter)
    if (qscript.baseQualityOutput != null) add(baseQualityCounter)
    if (qscript.readMappingQualityOutput != null) add(readQualityCounter)
    if (qscript.coverageOutput != null ) add(coverageCounter)
  \}
\}
