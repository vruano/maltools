package Maltools::Tools::Picard::MergeSamFiles;

use strict;
use warnings;
use Moose;
extends 'Maltools::DataTemplateTool';

our $INPUTS = { in_files => { type => 'bam', required => 1, multiple => 1 }, comments => {type => 'string', multiple => 1 } };
our $OUTPUTS = { out => { type => {name => 'bam', indexed => 1}, required => 1} };

sub interpreter {
  return '$SHELL';
}

1;

__DATA__
{ $input_arguments = join (" ", map { $_ eq "-" ? "INPUT=/dev/stdin" : "INPUT=$_" } @{$J->input("in_files")});
  $comment_arguments = join (" ", map { "COMMENT='$_'" } @{$J->input("comments")});
  $output_arguments = "OUTPUT=" . ($_ eq "-" ? '/dev/stdout': $J->output("out"));
  $generated_bai = $J->output("out"); $generated_bai =~ s/\.bam$//; $generated_bai .= ".bai";
  $output_bai = $J->output("out") . '.bai'; ''}
picard MergeSamFiles {$input_arguments} {$comment_arguments} {$output_arguments} CREATE_INDEX={$J->output("out") eq "-" ? 'false' : 'true'} MERGE_SEQUENCE_DICTIONARIES=true
{ $J->output("out") eq "-" ? '' : "mv $generated_bai $output_bai"  } 
