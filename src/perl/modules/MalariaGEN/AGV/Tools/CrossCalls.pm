package MalariaGEN::AGV::Tools::CrossCalls;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
   in => { type => 'vcf', mandatory => 1 },
   matrices => { type => 'file', mandatory => 1},
   report => { type => 'file', mandatory => 1},
   parents => { type => 'string', multiple => 1, mandatory => 1},
};

our $OUTPUTS = { 
   out => { type => 'vcf', mandatory => 1 },
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
  return 30 * 60;
}

sub interpreter {
   my ($self) = @_;
   return 'sh';
}

1;

__DATA__
{ 
  $in_vcf = $J->input('in');
  @parents  = @{$J->input('parents')};
  $matrices   = $J->input('matrices');
  $report  = $J->input('report');
  $out  = $J->output('out'); '' }

crossToVcf --in-vcf {$in_vcf} --matrix {$matrices} --report {$report} {join(" ",map { "-p $_" } @parents)} --output {$out}


