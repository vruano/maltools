package MalariaGEN::AGV::Tools::RecalAnalysis;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
   vcf => { type => { name => 'vcf' }, mandatory => 1 },
   before => { type => { name => 'bam', indexed => 1 }, mandatory => 1},
   after => { type => { name => 'bam', indexed => 1 }, mandatory => 1},
   ref => { type => 'fasta', mandatory => 1},
   baq => { type => 'bool', mandatory => 0, default => 1 },
   sample => { type => 'string', mandatory => 1},
};

our $OUTPUTS = { 
    out_dir => { type => 'file', mandatory => 1 },
};

sub job {
  my $self = shift;
  my $job = $self->SUPER::job(@_);
  $job->tmp_required(2000);
  return $job;
}

sub calculate_memory {
  return 2000;
}

sub calculate_cpu_ratio {
  return 2;
}

sub calculate_cpu_time {
  return 60 * 60;
}

sub interpreter {
   my ($self) = @_;
   return 'sh';
}

1;

__DATA__
{ 
  $vcf = $J->input('vcf');
  $ref  = $J->input('ref');
  $sample  = $J->input('sample');
  $before   = $J->input('before');
  $after  = $J->input('after');
  $baq  = $J->input('baq');
  $out  = $J->output('out_dir'); '' }
mkdir -p {$out}
bqErrorMatrix {$baq ? '--baq':'--nobaq'} --vcf {$vcf} --sample {$sample} --before {$before} --after {$after} --mat {$out}/ab.mat --bdo {$out}/b.dist --ado {$out}/a.dist --ref {$ref}


