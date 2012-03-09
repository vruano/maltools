package MalariaGEN::AGV::Tools::Bam::Coverage;

use Moose;
extends 'MalariaGEN::AGV::DataTemplateTool';

our $INPUTS = {
  in => { type => { name => 'bam', indexed => 1 }, mandatory => 1 },
  remove_duplicates => { type => 'bool', default => 0 },
  region_size => { type => 'num', default => 10000 },
};

our $OUTPUTS = {
  out => { type => "file" , mandatory => 1 },
};

sub interpreter {
  return '$SHELL';
}

__DATA__
{$in = $J->input('in');
 $out = $J->output('out');
 $region_size = $J->input('region_size');
 $rd = $J->input('remove_duplicates') ? '--remove-duplicates' : '';
 '' }
bam-coverage -i {$in} -o {$out} {$rd} -s {$region_size} || exit $?
