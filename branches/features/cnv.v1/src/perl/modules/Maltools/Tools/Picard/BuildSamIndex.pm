package Maltools::Tools::Picard::BuildSamIndex;

use strict;
use warnings;
use Moose;
extends 'Maltools::DataTemplateTool';

our $INPUTS = { in => { type => 'bam', required => 1 } };
                
our $OUTPUTS = { out => { type => 'file', required => 1} };

sub interpreter {
  return '$SHELL';
}

1;

__DATA__
{ $in = $J->input("in"); $in = '/dev/stdin' if $in eq '-'; 
  $out = $J->output("out"); $in = '/dev/stdout' if $in eq '-'; '' }
picard BuildBamIndex INPUT={$in} OUTPUT={$out}
