package MalariaGEN::AGV::Commands::DubiousRegions;

use Moose;
extends 'MalariaGEN::AGV::Command';

#has '+engine_name' => (default => 'local');

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'genotypes a set of samples';
}

sub hidden {
  return 1;
}

sub help_text {
   return "dubious-regions\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " dubious-regions [-i sample.bam ]+ -o output.intervals [--reference ref.fa]? "  
                         . "[--indels possible_indels.vcf | --indels-output detected_indels.vcf]?\n" .
          "\tYou can indicate a list of indels to include amongs the dubious regions with --indels\n" .
          "\totherwise the command will create its own list using GATK's SomaticIndelGenotyper\n" .
          "\tIf the latter you might as well recover the list of detected indels indicating and output\n" .
          "\tusing --indels_output.\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};

  my @inputs = ();
  my $output = "-";
  my $indels_input = undef;
  my $indels_output = undef;
  my $reference = undef;
  my @regions = ();
  my $class = "realigned";
  GetOptions("input|in|i=s@" => \@inputs, "output|out|o=s" => \$output,
             "reference|ref|r=s" => \$reference,
             "indels-input|indels=s" => \$indels_input,
             "indels-output=s" => \$indels_output);

  return $self->error_return("you must at least specify one sample") if $#inputs < 0;
  my @missing_inputs = grep { ! -f $_ } @inputs;
  return $self->error_return("the following input samples could not be reached or are not regular files: " . join(", ",@missing_inputs)) 
     unless $#missing_inputs < 0;
  return $self->error_return("you cannot specify and input and output indels files simultaneously")
    if defined($indels_input) && defined($indels_output);  
  
  my $rc = reference_config();
  
  $reference ||= $rc->file_name();
  my $indels_are_temp = !defined($indels_input) && !defined($indels_output);
  my $indels_file = $indels_input || $indels_output || $output . ".indels.vcf";

  my $engine = $self->resolve_engine;
  my $tool =  MalariaGEN::AGV::Tool->tool_by_name('gatk-dubiousRegions');
  my $job = $tool->job(
      inputs => { 
          ref => $reference, 
          samples => [@inputs] , 
          interval_list => $rc->file_name(extension => '.interval_list'),
#          (-f $indels_file) ? (indels_in => $indels_file) : (),
          scatter_count => 20,
      },
      outputs => { 
          regions_out => $output,
#          (! -f $indels_file) ? (indels_out => $indels_file) : (),
      }
  );
  unless ($engine->run_job($job)) {
    unlink $indels_file if $indels_are_temp;
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotyping with message: " . $job->error_message);
  }
  else {
    unlink $indels_file if $indels_are_temp;
    return $self->ok_return;
  }
}

1;
