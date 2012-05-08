package Maltools::Commands::CrossReportV2;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);
use Cwd qw(realpath);

use Maltools::Config qw(data_config 
  sanger_config reference_config jobs_config);

use Maltools::Engines::Sanger;
use Getopt::Long;
extends 'Maltools::Command';

sub hidden {
  return 1;
}

sub help_summary {
   return 'constructs cross accuracy report v2 from a set of cross matrices';
}

sub help_text {
   return "snps-filter\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " cross-report-v2 -i cross-matrices -c coverage-tables -o out-dir\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $manifest = undef;
  my $rc = reference_config();
  my $input = undef;
  my $output = undef;
  my $coverage = undef;
  my @parents= ();
  GetOptions("parent|p=s" => \@parents, "input|i=s" => \$input, "output|o=s" => \$output, "coverage|c=s" => \$coverage);

  $#parents == 1 or return $self->error_return("the number of parents provided must be two");
  -e $input && -d $input or return $self->error_return("you must provide an input cross matrix directory");
  -e $coverage && -d $coverage or return $self->error_return("the coverage table directory must exists and be a directory ($coverage)");
  -e $output && ! -d $output and return $self->error_return("the output directiory exists already but is not a directory");
  mkdir($output);
 
  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('crossReportV2');
  my $job = $tool->job(inputs => { in_dir => $input, coding_tables => $rc->coding_tables(version => 'v2'), 
               uniqueness_tables => $rc->uniqueness_tables, parents => \@parents,
               coverage_tables => $coverage } , outputs => { out_dir => $output });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
