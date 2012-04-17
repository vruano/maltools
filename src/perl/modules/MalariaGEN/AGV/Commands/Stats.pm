package MalariaGEN::AGV::Commands::Stats;

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
   return 'show statistics on a data object';
}


sub hidden {
   return 1;
}

sub help_text {
   my $self = shift;
   my $prog = $self->cl_name;
   my $summary = $self->help_summary;
   return <<EOT
Command:

  stats - $summary

Synopsis:

  $prog stats --bam input.bam -o output.txt

Description:

  Currently only alignemnt are supported and the statistics reported are those from samtools.
  However this shall extend to other data object such as vcf, sample, lanes etc...

EOT
;
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV = @{$params->{arguments}};
  my $input = undef;
  my $output = "-";
  GetOptionsFromArray(\@ARGV,"bam|i=s" => \$input,"o|output=s" => \$output);
  defined $input or return $self->error_return("you must provide an input alignment");
  -f $input or return $self->error_return("input alignment '$input' does not seem to exists or is not a regular file");  
  my $tool = MalariaGEN::AGV::Tool->tool_by_name("samtools-flagstat");
  my $job = $tool->job(inputs => { in => $input } , outputs => { out => $output });
  my $engine = $self->resolve_engine;
  $engine->run_job($job) 
	or $self->error_return("error executing the align command with exit code : " . $job->return_code,$job->error_message );
  return $self->ok_return;
}


1;
