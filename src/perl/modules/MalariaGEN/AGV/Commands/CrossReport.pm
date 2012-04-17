package MalariaGEN::AGV::Commands::CrossReport;

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

sub hidden {
  return 1;
}

sub help_summary {
   return 'constructs cross accuracy report material from a set of cross matrices';
}

sub help_text {
   return "snps-filter\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " cross-report -i cross-matrices -o out-dir\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $manifest = undef;
  my $rc = reference_config();
  my $input = undef;
  my $output = undef;
  GetOptions("input|i=s" => \$input, "output|o=s" => \$output);

  -e $input && -d $input or return $self->error_return("you must provide an input cross matrix directory");
  -e $output && ! -d $output and return $self->error_return("the output directiory exists already but is not a directory");
  mkdir($output);
 
  my $engine = $self->resolve_engine();
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('crossReport');
  my $job = $tool->job(inputs => { in_dir => $input, coding_tables => $rc->coding_tables, ref_uq => $rc->file_name(extension => '.uq') } , outputs => { out_dir => $output });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
