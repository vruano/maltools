package MalariaGEN::AGV::Tools::Samtools::Flagstat;

use strict;
use warnings;
use Moose;
extends 'MalariaGEN::AGV::DataTemplateTool';

our $INPUTS = { in => { type => 'bam', required => 1 } };
our $OUTPUTS = { out => { type => 'file', required => 1} };

sub interpreter {
  return '$SHELL';
}

sub calculate_cpu_time {
  return 2 * 60;
}


1;

__DATA__
samtools flagstat {$J->input("in")} > {$J->output("out")}

