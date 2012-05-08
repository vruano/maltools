package Maltools::Commands::CandidateSnps;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);

use Maltools::Config qw(data_config 
  sanger_config reference_config jobs_config);

use Getopt::Long;
extends 'Maltools::Command';

sub help_summary {
   return 'from vcf files calculate the candidate snplist';
}

sub hidden {
  return 1;
}

sub help_text {
   my $self = shift;
   my $program = $self->cl_name;
   my $summary = $self->help_summary;

   return <<EOM
Command:
  
  candidate-snps - $summary

Syntaxis:

  $program candidate-snps create -i input1.vcf -i input2.vcf ... -o output.snps [-r reference.fa]

    constructs a list of candidate snps from a set of Vcf files containing samtools denovo genotype calls.
    The default reference is used in non is provided.

  $program candidate-snps merge -i input1.snps -i input2.snps ... -o output.snps [-r reference.fa]

    constructs a list of candidate snps from a set of snps lists. It merges the content of the input list
    into a single list. The default reference is used if none is provided.

  $program candidate-snps props -i input.snps -o output.snps-props [-r reference.fa] [-a annotation.gff]
  
    constructs the snps property files used in some analyses downstream. The default reference and its default
    annotation file are used if none is provided.

EOM
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my @inputs = ();
  my $from_vcf = 0;
  my $reference = undef;
  my $annotation = undef;
  GetOptions(
    "input|in|i=s@" => \@inputs, 
    "output|out|o=s" => \$output,
    "reference|ref|r=s" => \$reference,
    "annotation|anno|a=s" => \$annotation);
  my $operation = shift(@ARGV) || "create"; 

  my @bad = grep { ! -e $_ } @inputs;
  return $self->error_return("some inputs are missing: " . join(", ",@bad)) if $#bad >= 0;

  @bad = grep { ! -f $_ } @inputs;
  return $self->error_return("some inputs are not regular files: " . join(", ",@bad)) if $#bad >= 0;

  @bad = grep { ! -r $_ } @inputs;
  return $self->error_return("some inputs cannot be read: ".join(", ",@bad)) if $#bad >= 0;
  
  return $self->error_return("you must specify and output") unless $output;
  return $self->error_return("the output if already existing must be a regular file") if -e $output && !-f $output;   

  my $config = Maltools::Config->new();

  $reference ||= $config->get_reference_sequences();
  $annotation ||= $config->get_reference_annotation();

  -f $reference or return $self->error_return("the reference file '$reference' cannot be open for reading");
  -f $annotation or return $self->error_return("the annotation file '$annotation' cannot be open for reading");

  my $engine = $self->resolve_engine();
  
  my $tool;
  my $outparam = "out";
  if ($operation eq "create") {
    $tool = Maltools::Tool->tool_by_name('vcfToVarLists');
    $outparam = "snpout";
  }
  elsif ($operation eq "merge") {
    $tool = Maltools::Tool->tool_b_name('mergeVarLists');
  }
  elsif ($operation eq "props" || $operation eq "properties") {
    $tool = Maltools::Tool->tool_by_name('varListsProperties');
  }
  my $job = $tool->job(inputs => { in=> \@inputs, ref => $reference, annotation => $annotation }, outputs => { $outparam => $output });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
