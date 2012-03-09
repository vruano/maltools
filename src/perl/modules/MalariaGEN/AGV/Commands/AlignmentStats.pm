package MalariaGEN::AGV::Commands::AlignmentStats;

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
   return "alignment-stats\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " alignment-stats -i input -o output\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $input = "-"; 
  my $output = "-";
  GetOptions("i|input=s" => \$input,"o|output=s" => \$output);
  -f $input or return $self->error_return("input alignment '$input' does not seem to exists or is not a regular file");  
  my $tool = MalariaGEN::AGV::Tool->tool_by_name("samtools-flagstat");
  my $job = $tool->job(inputs => { in => $input } , outputs => { out => $output });
  my $engine = $self->resolve_engine;
  $engine->run_job($job) 
	or $self->error_return("error executing the align command with exit code : " . $job->return_code,$job->error_message );
  return $self->ok_return;
}


1;
