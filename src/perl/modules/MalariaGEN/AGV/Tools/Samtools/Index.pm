package MalariaGEN::AGV::Tools::Samtools::Index;

use strict;
use warnings;
use Moose;
extends 'MalariaGEN::AGV::DataTemplateTool';

our $INPUTS = { in => { type => 'bam', required => 1, indexed => 0 } };
                
our $OUTPUTS = { out => { type => 'file', required => 1} };

sub interpreter {
  return '$SHELL';
}

1;

__DATA__
{ $in = $J->input("in"); $in = '/dev/stdin' if $in eq '-'; 
  $out = $J->output("out"); $in = '/dev/stdout' if $in eq '-'; '' }
samtools index {$in}
{ $out ne ($in . ".bai") ? "mv ${in}.bai $out" : "" }
