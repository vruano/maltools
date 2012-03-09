package MalariaGEN::AGV::Commands::Merge;

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
   my $pgv = $_[0]->cl_name;
   my $summary = $_[0]->help_summary;

   return <<EOM
Summary:
 
  merge - $summary

Synapsys:

  $pgv merge -i input1.bam -i input2.bam ... -o output.bam

EOM
;
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my @inputs = ();
  my $output = "-";
  GetOptions("i|input=s@" => \@inputs,"o|output=s" => \$output);
  push @inputs,@ARGV;
  
  my $tool = MalariaGEN::AGV::Tool->tool_by_name("picard-mergeSamFiles");
  my $rc = reference_config(); 
  my @real_inputs = map { -e $_ ? realpath($_) : $_ } @inputs;
  my $job = $tool->job(inputs => { ref => $rc->file_name, in_files => \@real_inputs } , outputs => { out => $output });
  my $engine = $self->resolve_engine;
  $engine->run_job($job) 
	or $self->error_return("error executing the align command with exit code : " . $job->return_code,$job->error_message );
  return $self->ok_return;
}


1;
