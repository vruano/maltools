package MalariaGEN::AGV::Commands::MarkDuplicates;

use strict;
use warnings;
use File::Basename qw(dirname);
use Getopt::Long;

use Moose;

extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'mark duplicate read in a sam/bam file';
}

sub hidden {
  return 1;
}

sub help_text {
   return "map\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " mark-duplicates -i input-bam -o output-bam\n";
}

sub execute {
  my ($self,$site,$params) = @_;

  local @ARGV = @{$params->{arguments}};
  my $input = "-";
  my $output = "-";
  my $metrics = "";
  my $remove = 0;
  GetOptions("input|i=s"=> \$input,"output|o=s" => \$output, "metrics|m=s" => \$metrics, "remove|R!" => \$remove);

  return $self->error_return("you must specify a existing regular file as input '$input'") unless -f $input;
  
  my $out_dir = $output eq "-" ? "" : dirname($output);
  unless ($output eq "-" || -d $out_dir) {
    return $self->error_return("the output file containing directory does not exist or is not in fact a directory");
  }
  my $engine = $self->resolve_engine;
  my $tool = MalariaGEN::AGV::Tool->tool_by_name("picard-markDuplicates");
  my $job = $tool->job(inputs => { in => $input, remove => $remove } , outputs => { out => $output, ($metrics ? (metrics => $metrics) : ()) });
  $engine->run_job($job) 
    or return $self->error_return($job->error_message);
  return $self->ok_return;
}

1;
