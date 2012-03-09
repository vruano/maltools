package MalariaGEN::AGV::Tools::Picard::SortSam;

use strict;
use warnings;
use Moose;
extends 'MalariaGEN::AGV::DataTemplateTool';

our $INPUTS = { in => { type => 'bam', required => 1 },
                order => { type => 'string', default => 'coordinate' }};
our $OUTPUTS = { out => { type => { name => 'bam', indexed => 1}, required => 1} };

sub interpreter {
  return '$SHELL';
}

1;

__DATA__
{$in = $J->input("in"); $in = "/dev/stdin" if $in eq "-";
 $out = $J->output("out"); $out = "/dev/stdout" if $out eq "-";
 $out_bai = $out; $out_bai =~ s/\.bam/.bai/; '' }
picard SortSam INPUT={$in} OUTPUT={$out} CREATE_INDEX={$out eq "/dev/stdout"?'false':'true'} SORT_ORDER={$J->input("order")} || exit $?
{  ($out ne "/dev/stdout") ? "mv $out_bai ${out}.bai" : '' }
