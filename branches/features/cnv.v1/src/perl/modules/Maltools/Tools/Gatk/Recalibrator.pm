package Maltools::Tools::Gatk::Recalibrator;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
    ref => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
    in => { type => {name => 'bam', indexed => 1}, mandatory => 1 },
    mask => { type => { name => 'file' }, mandatory => 1 },
    mask_format => { type => 'string', mandatory => 1 },
    recal_in => { type => { name => 'bam', indexed => 1}, default => undef},
};

our $OUTPUTS = { 
    out => { type => { name => 'bam', indexed => 1 }, mandatory => 1 },
    report => { type => { name => 'file' }, mandatory => 0 },
    recal_file => { type => 'file' , mandatory => 0 },
};


sub calculate_cpu_ratio {
  return 1;
}

sub calculate_memory {
  return 2000;
}

sub calculate_cpu_time {
  return 2 * 60 * 60;
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;

__DATA__

{ 
  $ref = $J->input("ref");
  $in_bam = $J->input("in");
  $recal_bam = $J->input("recal_in") || $J->input("in");
  $out_bam = $J->output("out");
  $out_rep = $J->output("report");
  $out_recal = $J->output("recal_file");
  $rep_tmp = $J->tempdir();
  $mask = $J->input("mask");
  $mask_format = $J->input("mask_format"); '' }

mkdir -p {$rep_tmp}/before || exit 1;
mkdir -p {$rep_tmp}/after || exit 2;

gatk --memory {$J->memory - 100} -T CountCovariates --default_platform ILLUMINA -standard -R {$ref} -I {$recal_bam} -recalFile {$rep_tmp}/before/output.recal_data.csv -B:mask,{$mask_format} {$mask} || exit 3
gatk --memory {$J->memory - 100} -T TableRecalibration --default_platform ILLUMINA -R {$ref} -I {$in_bam} -recalFile {$rep_tmp}/before/output.recal_data.csv --out {$out_bam} || exit 4
picard BuildBamIndex INPUT={$out_bam} OUTPUT={$out_bam}.bai || exit 5

{ $out_recal ? "cp ${rep_tmp}/before/output.recal_data.csv $out_recal" : "" }

{ $out_rep ? '' : "rm -Rf ${rep_tmp}; exit 0" }

gatk --memory {$J->memory - 100} -T CountCovariates --default_platform ILLUMINA -standard -R {$ref} -I {$out_bam} -recalFile {$rep_tmp}/after/output.recal_data.csv -B:mask,{$mask_format} {$mask} || exit 6

gatk-ac --memory {$J->memory - 100} -recalFile {$rep_tmp}/before/output.recal_data.csv -outputDir {$rep_tmp}/before || exit 7
gatk-ac --memory {$J->memory - 100} -recalFile {$rep_tmp}/after/output.recal_data.csv -outputDir {$rep_tmp}/after || exit 8

pushd {$rep_tmp} || exit 9
tar czf {$out_rep} . || exit 10
popd || exit 11

exit 0
