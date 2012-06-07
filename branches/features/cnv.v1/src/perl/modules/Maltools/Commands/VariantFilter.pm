package Maltools::Commands::VariantFilter;

use strict;
use warnings;
use Moose;

use Maltools::Config qw(data_config 
       sanger_config reference_config jobs_config);

use Maltools::Engines::Sanger;
use Getopt::Long;
extends 'Maltools::Command';

sub help_summary {
   return 'filters variants from a vcf file';
}


sub hidden {
   return 1;
}

sub help_text {
   return "genotype\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " variant-filter FILTER1 FILTER2 FILTER3 ... [-i input-file] [-o output-file]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my $input = "-";
  GetOptions("i=s" => \$input, "o=s" => \$output);
   
  my $dc = data_config();
  my $rc = reference_config();
  my $tool = Maltools::Tool->tool_by_name('vcf-filter');
  my $engine = Maltools::Engines::Sanger->new();
  my $job = $tool->job(inputs => { in => $input, filters => [@ARGV], op => 'add' }, outputs => { out => $output});
  $engine->run_job($job);
  return $self->ok_return;
}


1;
