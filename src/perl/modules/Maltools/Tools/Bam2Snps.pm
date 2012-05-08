package Maltools::Tools::Bam2Snps;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
   "candidates" => { type => 'file', mandatory => 0},
   "minq" => { type => 'num', mandatory => 0, default => 27},
   "ref" => { type => 'fasta', mandatory => 1},
   "in" => { type => { name => 'bam', indexed => 1 },  mandatory => 1},
};

our $OUTPUTS = { 
    out => { type => { name => 'file' }, mandatory => 1 },
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 1000;
}

sub calculate_cpu_time {
  return 10 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;

__DATA__
{ 
  $snplist_argument = $J->input('candidates') ? '--snplist=' . $J->input('candidates') : '';
  $ref  = $J->input('ref');
  $in   = $J->input('in');
  $minq = $J->input('minq');
  $out  = $J->output('out'); '' }
samtools mpileup -Bf {$ref} {$in} | bam2snps --mode=mpileup {$snplist_argument} --snpout={$out} --minq={$minq}

