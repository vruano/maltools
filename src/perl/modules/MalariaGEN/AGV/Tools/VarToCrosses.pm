package MalariaGEN::AGV::Tools::VarToCrosses;

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
   in => { type => 'vcf', mandatory => 1 },
   parents => { type => 'string', multiple => 1, mandatory => 1},
   invariants => { type => 'string', default => 'no' },
   includes => { type => 'string', multiple => 1, mandatory => 1},
   excludes => { type => 'string', multiple => 1, mandatory => 1}
};

our $OUTPUTS = { 
   "out_dir" => { type => 'file' }
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 100;
}

sub calculate_cpu_time {
  return 5 * 60;
}

sub interpreter {
  my ($self) = @_;
  return '$SHELL';
}

1;

__DATA__
{ $in  = $J->input('in');
  $out = $J->output('out_dir');
  $exclude_invariants = $J->input('invariants');
  $ei_params = $exclude_invariants eq 'no' ? "" : "--excl-invs=$exclude_invariants";
  $includes = $J->input('includes');
  $excludes = $J->input('excludes');
  $parents = $J->input('parents'); '' }
echo $PATH
`which vcfToCrosses`
vcfToCrosses -i {$in} -o {$out} {join(" ", map { "-p $_" } @$parents)} {join(" ", map { "--include $_" } @$includes)} {join(" ",map { "--exclude $_" } @$excludes)} {$ei_params}

