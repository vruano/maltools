package MalariaGEN::AGV::Tools::Bam2Cov;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
   ref => { type => 'fasta', mandatory => 1},
   in => { type => { name => 'bam', indexed => 1 },  mandatory => 1},
};

our $OUTPUTS = { 
    out => { type => { name => 'file' }, mandatory => 1 },
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 100;
}

sub calculate_cpu_time {
  return 30 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;

__DATA__
{ 
  $ref  = $J->input('ref');
  $in   = $J->input('in');
  $out  = $J->output('out'); '' }
bam2cov --reference {$ref} {$in} > {$out}
