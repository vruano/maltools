package MalariaGEN::AGV::Commands::Setup;

use base 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'setup the AGV pipeline initializing reference and database entries if not existent';
}

sub help_text {
   return "setup:\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " setup\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  my $rc = reference_config();
  my $bwa_index = $rc->file_name(extension => ".fa.sa");
  -e $bwa_index or $self->_bwa_index($rc);
  return $self->ok_return;
}


sub _bwa_index {
  my ($self,$rc) = @_;
  my $job = MalariaGEN::AGV::Tool->tool_by_name("bwa-index")->job(inputs => {
         in => $rc->file_name(extension => ".fa")});
  my $engine = jobs_config()->engine;
  $engine->run_job(job => $job);
  return $self->ok_return;
}



1;
