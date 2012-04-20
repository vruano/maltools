package MalariaGEN::AGV::Tools::CrossReport;

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
   in_dir => { type => 'file', mandatory => 1 },
   ref_uq => { type => 'file', mandatory => 1 },
   coding_tables => { type => 'file', mandatory => 1},
   parent0 => { type => 'string', default => 'parentA' },
   parent1 => { type => 'string', default => 'parentB' },
};

our $OUTPUTS = { 
   "out_dir" => { type => 'file', mandatory => 1 }
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 4000;
}

sub calculate_cpu_time {
  return 60 * 60;
}

sub scripts_dir {
   my ($self) = @_;
   return realpath(catfile($ENV{MALTOOLS_HOME},'resources','cross_report','Scripts'));
}

sub interpreter {
  my ($self) = @_;
  return 'python';
}

1;

__DATA__
{ $scripts_dir = $T->scripts_dir;
  $parent0 = $J->input('parent0');
  $parent1 = $J->input('parent1');
  $out_dir = $J->output('out_dir');
  $ref_uq = $J->input('ref_uq');
  $coding_tables = $J->input('coding_tables');
  $in_dir = $J->input('in_dir'); '' }
import os
import readline
import string
import sys
import random
from os import system
import re

## Set the key variables:
recomb_rate = 0.0001;
geno_rate   = 0.001;

## Compile the C++ code:
#line = "g++ $(pkg-config --cflags gsl) -c hmm.cpp";
#os.system(line);

#line = "g++ hmm.o $(pkg-config --libs gsl) -o hmm"
#os.system(line);

os.system("cp -Rf {$in_dir}/* {$out_dir}");

## Run the HMM to produce Viterbi (color) paths:
line = "hmm {$out_dir}/ 0.00001 0.0001";
os.system(line);

os.system("mkdir -p {$out_dir}/Output");
os.system("mkdir -p {$out_dir}/Figures");
## R script sequence to create color figures of cross results

#>	parse out figure path


fig_path     = "{$out_dir}/Figures"

cross_1      = "{$parent0}"
cross_2      = "{$parent1}"


#>	execute script

command = "cat {$scripts_dir}/calc.errors.r | R --no-save --args {$out_dir}/ {$out_dir}/Figures/ {$parent0} {$parent1}";
print(command);
os.system(command)

## R script to create pipeline data files

#>	parse out output path

out_path     = "{$out_dir}/Output";

#>	execute script

command = "cat {$scripts_dir}/plot.cnc.uniq.r | R --no-save --args {$out_dir}/ {$out_dir}/Output/ {$ref_uq} {$coding_tables}";
print(command)
os.system(command);

