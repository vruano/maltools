package MalariaGEN::AGV::Tools::Bwa2Som;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
   "snps" => { type => 'file', mandatory => 1}, 
   "ref" => { type => 'fasta', mandatory => 1},
   "in" => { type => { name => 'bam', indexed => 1 },  mandatory => 1},
};

our $OUTPUTS = { 
    out => { type => { name => 'bam', indexed => 1 }, mandatory => 1 },
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub job {
  my $self = shift;
  my $job = $self->SUPER::job(@_);
  $job->tmp_required(2000);
  return $job;
}

sub calculate_memory {
  return 2000;
}

sub calculate_cpu_time {
  return 45 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return "perl";
}

1;

__DATA__
{ 
  $snps = $J->input('snps');
  $ref  = $J->input('ref');
  $in   = $J->input('in');
  $out  = $J->output('out'); '' }

use Getopt::Long qw(:config ignore_case);
use File::Temp qw(tempfile);
use File::Spec::Functions qw(catfile);
use Cwd qw(cwd);
use strict;
use warnings;

my $snpList = "{ $snps || 'all' }";
my $reference = "{ $ref }";
my $input = "{ $in }";
my $output = "{ $out }";

my $bwa2som_arguments = "--bam=$input --ref=$reference --indels=all --out=$output --snps=$snpList --fix";

`bwa2som $bwa2som_arguments 2> /dev/stderr`;
if ($?) \{
  die "failed in exectuing the bwa2som step";
\}

`samtools index $output`;
if ($?) \{
  die "failed to create an index of the alignment '$output'";
\}
