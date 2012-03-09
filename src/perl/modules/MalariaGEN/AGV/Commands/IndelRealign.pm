package MalariaGEN::AGV::Commands::IndelRealign;

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
   return 'realigns sample BAM files around indels';
}

sub help_text {
   return "indel-realign\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " indel-realign [-i input-bam | < input-bam ] --intervals intervals-file [-o output-file | > output-file]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my $input = "-";
  my $intervals = undef;
  GetOptions("i=s" => \$input, "o=s" => \$output, "intervals=s" => \$intervals);
  $intervals or return $self->error_return("you must provide an internvals file");
  -e $intervals or return $self->error_return("cannot access '$intervals'");
  -f $intervals or return $self->error_return("'$intervals' is not a regular file"); 
  my $engine = $self->resolve_engine();
  my $ref = reference_config->file_name();
  my $real_input = $input;
  if ($input eq "-") {
     $real_input = $engine->tempfile(); 
     print STDERR "Coping input into a file '$real_input' ...";
     copy(\*STDIN,$real_input) or return $self->error_return("could not copy input in temporal file with name '$real_input'");
  } 
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-indelRealigner');
  my $job = $tool->job(inputs => { ref => $ref, inputs => [$real_input], intervals => $intervals }, outputs => { out => $output});
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
  
}


1;
