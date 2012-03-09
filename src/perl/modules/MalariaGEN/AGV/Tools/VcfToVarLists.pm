package MalariaGEN::AGV::Tools::VcfToVarLists;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use vars qw(%ENV);
use File::Basename qw(dirname);
use Text::Template;
use IO::File;
use Cwd qw(realpath);
use JSON::XS;
use File::Spec::Functions qw(catfile file_name_is_absolute);
use POSIX;

our $INPUTS = {
   "in" => { type => 'vcf', multiple => 1, mandatory => 1},
};


our $OUTPUTS = { 
   "snpout" => { type => 'file', mandatory => 1 },
   "indelout" => { type => 'file', mandatory => 0 },
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 2000;
}

sub calculate_cpu_time {
  return 4 * 60 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}


1;

__DATA__
{  
 $snpout = $J->output('snpout') || '/dev/null';
 $indelout = $J->output('indelout') || '/dev/null';
 @inputs = @{$J->input('in')};
 $inputs_string = join(" ",map { "-i $_" } @inputs);
 '' }

vcfToVarLists {$inputs_string} --snpout {$snpout} --indelout {$indelout}


