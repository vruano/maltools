package Maltools::Tools::Gatk::IndelRealigner;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
    ref => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
    inputs => { type => { name => 'bam', indexed => 1 }, multiple => 1, mandatory => 1 },
    intervals => { type => 'intervals', multiple => 0, mandatory => 1 },
};

our $OUTPUTS = { 
    out => { type => { name => 'bam', indexed => 1 }, mandatory => 1 }
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 6000;
}

sub calculate_cpu_time {
  return 60 * 60 * 3;
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;

__DATA__
{ $mem = $J->memory;
  $cpu_params = $J->cpu_count == 1 ? "":"-nt " . $J->cpu_count . " ";
  $ref = $J->input("ref");
  $inputs_params =join(" ",map {"-I $_"} @{$J->input("inputs")} );
  $intervals = $J->input("intervals");
  $out = $J->output("out"); 
  $output_params = $out eq "-" ? '' : "-o $out" ; '' }
gatk --memory {$mem} {$cpu_params} -T IndelRealigner -R {$ref} {$inputs_params} -targetIntervals {$intervals} {$output_params}
