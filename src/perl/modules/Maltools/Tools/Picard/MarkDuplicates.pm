package Maltools::Tools::Picard::MarkDuplicates;

use strict;
use warnings;
use Moose;
extends 'Maltools::DataTemplateTool';

our $INPUTS = { in => { type => 'bam', required => 1 },
                remove => { type => 'bool', default => 0 } };
our $OUTPUTS = { out => { type => 'bam', indexed => 1, required => 1 }, metrics => { type => 'file', required => 0 } };

sub interpreter {
  return '$SHELL';
}

1;

__DATA__
{ $in = $J->input("in"); $in = '/dev/stdin' if $in eq "-";
  $metrics_out = $J->output("metrics") || '/dev/null';
  $remove = $J->input("remove") ? "true" : "false";
  $out = $J->output("out"); $out = '/dev/stdout' if $out eq "-"; '' }
picard MarkDuplicates INPUT={$in} OUTPUT={$out} REMOVE_DUPLICATES={$remove} METRICS_FILE={$metrics_out} AS=true CREATE_INDEX=true
