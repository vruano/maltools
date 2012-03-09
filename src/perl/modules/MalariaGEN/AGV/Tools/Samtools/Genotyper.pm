package MalariaGEN::AGV::Tools::Samtools::Genotyper;

use strict;
use warnings;
use Moose;
extends 'MalariaGEN::AGV::DataTemplateTool';

our $INPUTS = { in => { type => 'bam' , multiple => 1, mandatory => 1 }, 
                reference => { type => { name => 'fasta', indexed => 1} , mandatory => 1},
                max_depth => { type => 'num' , default => 100 }}; 
our $OUTPUTS = { out => { type => 'vcf', mandatory => 1} };

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
  $tmp_bcf = $J->tempfile("vcf_XXXX");
  $tmp_mpileup = $J->tempfile("mpu_XXXX");
  $out = $J->output('out');
  @samples = @{$J->input('in')};
  $max_depth = $J->input('max_depth');
'' }
samtools mpileup -C50 -uf {$ref} {join(" ",@samples)} > {$tmp_mpileup}  || exit 1;
bcftools view -bvcg {$tmp_mpileup} > {$tmp_bcf} || exit 2;
rm {$tmp_mpileup} || exit 4;
(bcftools view {$tmp_bcf} | vcfutils.pl varFilter -D{$max_depth} > {$out}) || exit 3;
rm {$tmp_bcf}

