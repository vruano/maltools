package MalariaGEN::AGV::Commands::IndelDiscovery;

use strict;
use warnings;
use Moose;

use File::Copy qw(copy);

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);

use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long;
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'list regions that might contain indels';
}

sub hidden {
  return 1;
}

sub help_text {
   return "indel-discovery\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " indel-discovery [ -i input.bam ]*  [-o output-file ]? \n" .
          "\t\tIf no input bam file is explicitly provided it would take the standard input instead\n" .
          "\t\tand similarly with the standard output\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my $input = "-";
  GetOptions("i=s" => \$input, "o=s" => \$output);
  my $engine = $self->resolve_engine();
  my $ref = reference_config->file_name();
  my $real_input = $input;
  if ($input eq "-") {
     $real_input = $engine->tempfile(); 
     print STDERR "Coping input into a file '$real_input' ...";
     copy(\*STDIN,$real_input) or return $self->error_return("could not copy input in temporal file with name '$real_input'");
  } 
  else {
    -e $input or return $self->error_return("cannot reach input file '$input'");
    -f $input or return $self->error_return("input file '$input' is not a regular file");
  }
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-realignerTargetCreator');
  my $job = $tool->job(inputs => { ref => $ref, inputs => [$real_input]}, outputs => { out => $output});
  my $result = $engine->run_job($job);
#  print STDERR "Result $result " . $job->return_code . "\n";
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
