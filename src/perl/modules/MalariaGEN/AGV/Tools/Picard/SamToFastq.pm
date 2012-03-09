package MalariaGEN::AGV::Tools::Picard::SamToFastq;

use strict;
use warnings;
use Moose;
extends 'MalariaGEN::AGV::DataTemplateTool';

our $INPUTS = { in => { type => 'bam', required => 1 } };
our $OUTPUTS = { forward => { type => 'fastq', required => 1}, reverse => { type => 'fastq', required => 1} };

sub interpreter {
  return '$SHELL';
}

1;

__DATA__
picard SamToFastq INPUT={$J->input("in")} FASTQ={$J->output("forward")} SECOND_END_FASTQ={$J->output("reverse")}
