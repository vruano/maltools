package MalariaGEN::AGV::Commands::SortAlignment;

use base 'MalariaGEN::AGV::Command';

use Cwd qw(realpath);
use strict;
use warnings;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'merge read alignments from several sam/bam files into one';
}

sub help_text {
   return "sort-alignment\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " sort-alignment [-i input | < input] [-o output | > output]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $input = "-"; 
  my $output = "-";
  GetOptions("i|input=s" => \$input,"o|output=s" => \$output);
  
  my $tool = MalariaGEN::AGV::Tool->tool_by_name("picard-sortSam");
  my $job = $tool->job(inputs => { in => $input } , outputs => { out => $output });
  my $engine = $self->resolve_engine;
  $engine->run_job($job) 
	or $self->error_return("error executing the align command with exit code : " . $job->return_code,$job->error_message );
  return $self->ok_return;
}


1;
