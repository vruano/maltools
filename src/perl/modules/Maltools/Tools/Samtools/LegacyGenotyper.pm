package Maltools::Tools::Samtools::LegacyGenotyper;

use strict;
use warnings;
use Moose;
extends 'Maltools::DataTemplateTool';

our $INPUTS = { in => { type => 'bam' , multiple => 1, mandatory => 1 }, 
                reference => { type => { name => 'fasta', indexed => 1} , mandatory => 1},
                max_depth => { type => 'num' , default => 100 }}; 
our $OUTPUTS = { out => { type => 'file', mandatory => 1} };

sub interpreter {
  return '$SHELL';
}

sub calculate_cpu_time {
  return 2 * 60 * 60;
}

sub calculate_memory {
  return 2000;
}


1;

__DATA__
{
  $ref = $J->input('reference');
  $tmp_pileup = $J->tempfile("pu_XXXX");
  $out = $J->output('out');
  @samples = @{$J->input('in')};
  $max_depth = $J->input('max_depth');
'' }
samtools pileup -vcf {$ref} {join(" ",@samples)} > {$tmp_pileup}  || exit 1;
samtools.pl varFilter {$tmp_pileup} | awk '$6>=20' | awk '(index($3,"*") == 0)\{ print $1"\t"$2"\t"$3"\t"$4 \}'  > {$out} 

