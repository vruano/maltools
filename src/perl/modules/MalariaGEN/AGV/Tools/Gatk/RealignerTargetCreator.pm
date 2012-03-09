package MalariaGEN::AGV::Tools::Gatk::RealignerTargetCreator;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

  
our $INPUTS = {
    ref => {  type => 'fasta', indexed => 1,  mandatory => 1 },
    inputs => { type => {name => 'bam', indexed => 1}, multiple => 1, mandatory => 1 },
};

our $OUTPUTS = { 
    out => { type => 'intervals', mandatory => 1 }
};


sub calculate_cpu_ratio {
  return 1;
}

sub calculate_memory {
  return 1000;
}

sub calculate_cpu_time {
  return 20000;
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;

__DATA__
{
 $mem = $J->memory;
 $cpu_params =  $J->cpu_count == 1 ? "":"-nt " . $J->cpu_count . " ";
 $ref = $J->input("ref");
 $input_params = join(" ",map {"-I $_"} @{$J->input("inputs")} ); 
 $out = $J->output("out");
 $out_params = $out eq "-" ? "" : "-o $out"; "" }
echo gatk --memory {$mem} {$cpu_params} -T RealignerTargetCreator -R {$ref} {$input_params} {$out_params} > /dev/stderr
gatk --memory {$mem} {$cpu_params} -T RealignerTargetCreator -R {$ref} {$input_params} {$out_params}

