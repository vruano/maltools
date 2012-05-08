package Maltools::Tools::Gatk::DubiousRegions;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) }
  
our $INPUTS = {
    ref => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
    samples => { type => { name => 'bam', indexed => 1}, multiple => 1, mandatory => 1 },
    interval_list => { type => 'file', mandatory => 1 },
#    indels_in => { type => 'file', mandatory => 0 },
    scatter_count => { type => 'num' , default => 1}, 
};

our $OUTPUTS = { 
#    indels_out => { type => 'vcf', mandatory => 0 },
    regions_out => { type => 'intervals', mandatory => 1 },
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
  return 1000;
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
     . '-jobQueue long '
     . '-bsub -run '
#     . '-indels {$J->output("indels_out")} ' 
     . '-regions {$J->output("regions_out")} ' );
}


1;


__DATA__
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An example building on the intro ExampleCountReads.scala.
 * Runs an INCOMPLETE version of the UnifiedGenotyper with VariantEval and optional VariantFiltration.
 */
class DubiousRegions extends QScript \{
  
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R", required=false)
  var referenceFile: File = new File("{$J->input("ref")}")
 // { scalar(@{$J->input("samples")}) }
  @Input(doc="Bam files to inspect for dubious regions.", shortName="I",required=false)
  var bamFiles: List[File] = List({  join(", ",map { "new File(\"" . $_ . "\")" }  @{$J->input("samples")} ) })

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = new File("{ $J->input("interval_list") }")

//  @Output(fullName="indelsOutput", doc="output indel genotyped file",shortName="indels", required=true)
//  var indelsOutput: File = _
 
  @Output(fullName="regionsOutput", doc="output dubious regions file",shortName="regions", required=true)
  var regionsOutput : File = _
  
  @Argument(doc="Memory to be allocated for each process running this script, expressed in MB", shortName="M", required=false)
  var memory: Int = 1000

  @Argument(doc="Number of processes to scatter the work upon", shortName="P",required=false)
  var procs: Int = {$J->input("scatter_count")} 

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait GeneralArguments extends CommandLineGATK \{
    this.reference_sequence = qscript.referenceFile
    this.intervals = List(qscript.intervals)
    this.input_file = qscript.bamFiles
    this.memoryLimit = qscript.memory
  \}

  def script = \{
    // Create the four function that we can run.
//    val indelGenotyper = new SomaticIndelDetector with GeneralArguments
    val realignerTargetCreator = new RealignerTargetCreator with GeneralArguments
//    indelGenotyper.scatterCount = qscript.procs
//    indelGenotyper.out = qscript.indelsOutput
//    indelGenotyper.unpaired = true
//    realignerTargetCreator.rodBind :+= RodBind("indels", "VCF", indelGenotyper.out);
    realignerTargetCreator.out = qscript.regionsOutput
    realignerTargetCreator.scatterCount = qscript.procs
//    add(indelGenotyper,realignerTargetCreator)
    add(realignerTargetCreator)
  \}
\}
