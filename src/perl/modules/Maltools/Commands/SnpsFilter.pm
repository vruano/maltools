package Maltools::Commands::SnpsFilter;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);
use Cwd qw(realpath);

use Maltools::Config qw(data_config 
  sanger_config reference_config jobs_config);

use Maltools::Engines::Sanger;
use Getopt::Long;
extends 'Maltools::Command';

sub help_summary {
   return 'from snps list produces VCF with the appropriate Variant filters';
}

sub hidden {
   return 1;
}

sub help_text {
   return "snps-filter\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " snps-filter --manifest manifest-file\n"; 
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $manifest = undef;
  my $calculate_cutoffs = 0;
  my $rc = reference_config();
  GetOptions("manifest|man=s" => \$manifest,"cutoffs!" => \$calculate_cutoffs);
   
  $manifest or return $self->error_return("you must provide a manifest file");
  $manifest = realpath($manifest);
  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name($calculate_cutoffs ? 'covCutoffs' : 'varFilter');
  my $job = $tool->job(inputs => { manifest => $manifest });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
