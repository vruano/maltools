package MalariaGEN::AGV::Commands::CombineVariants;

use strict;
use warnings;
use Moose;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long qw(:config pass_through no_getopt_compat no_ignore_case);
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'combines diverse variant VCF files into one';
}

sub help_text {
   return "combine-variants\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " combine-variants -i input1.vcf -i input2.vcf ... -o output.vcf \n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my @inputs = ();
  GetOptions("input|i=s@" => \@inputs, "o=s" => \$output);
  $output ne "-" or return $self->error_return("you must specify a file for output");
  $#inputs >= 0 or return $self->error_return("you must specify at least one variant file");
  grep { ! -f $_ } @inputs and 
     return $self->error_return("there seems to be some input files inexistent or that are not regular files: " . join(", ", grep { ! -f $_ } @inputs));
  
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-combineVariants');
  my $engine = $self->resolve_engine();
  my $job = $tool->job(inputs => { in => \@inputs }, outputs => { out => $output});
  my $result = $engine->run_job($job);
  return $self->ok_return if $result;
  return $self->error_return("error combining variant files with exit code " . $job->return_code,$job->error_message());
}


1;
