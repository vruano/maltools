package MalariaGEN::AGV::Commands::Cp;

use base 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'copy files from one data location to another';
}

sub help_text {
   return "cp\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " cp org dest\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $org = shift @ARGV or die "you need to specify origin and destination\n";
  my $dest = shift  @ARGV or die "you need to specify destination\n";
  my $copy = MalariaGEN::AGV::Tool->tool_by_name("copy");
  my $job = $copy->job(inputs => { src => $org },
                       outputs => { dest => $dest });
  my $jc = jobs_config();
  my $engine = $jc->engine();
  $engine->run_job($job);
  return $self->ok_return;
}


1;
