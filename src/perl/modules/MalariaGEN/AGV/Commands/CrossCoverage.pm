package MalariaGEN::AGV::Commands::CrossCoverage;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);
use Cwd qw(realpath);

use MalariaGEN::AGV::Config qw(data_config 
  sanger_config reference_config jobs_config);

use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long;
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'constructs coverage tables digested to be used by cross analysis';
}

sub help_text {
   return "snps-filter\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " cross-coverage -i coverage-file1 -i coverage-file2 -i coverage-file3 ... -o out-dir\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $manifest = undef;
  my $rc = reference_config();
  my @input = ();
  my $output = undef;
  GetOptions("input|i=s@" => \@input, "output|o=s" => \$output);

  $#input >= 0 or die "you must provide at least one coverage input file";
  my @missing_inputs = grep {! -e $_ } @input;
  $#missing_inputs < 0 or die "some coverage input files are missing: " . join(", ", @missing_inputs );

  -e $output && ! -d $output and return $self->error_return("the output directiory exists already but is not a directory");
  mkdir($output);
 
  my $engine = $self->resolve_engine();
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('crossCoverage');
  my $job = $tool->job(inputs => { in_files => \@input } , outputs => { out_dir => $output });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
