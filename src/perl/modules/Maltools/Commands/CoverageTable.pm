package Maltools::Commands::CoverageTable;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);

use Maltools::Config qw(data_config 
  sanger_config reference_config jobs_config);

use Maltools::Engines::Sanger;
use Getopt::Long;
extends 'Maltools::Command';

sub hidden {
  return 1;
}

sub help_summary {
   return 'filters an alignment reads using the Snp-o-Matic approach (Bwa2Som)';
}

sub help_text {
   return "coverage-table\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " coverage-table -i input-bam -o output.tsv [--reference reference]?\n"; 
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $input = undef;
  my $reference = undef;
  my $rc = reference_config();
  GetOptions("input|in|i=s" => \$input, "output|out|o=s" => \$output, "reference|r=s" => \$reference);
  $reference ||= $rc->file_name();
   
 
  -f $input or return $self->error_return("the input alignment '$input' is unreachable or is not a regular file");
  -f $reference or return $self->error_return("the reference '$reference' is unreachable or is not a regular file");
   
  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('bam2Cov');
  my $job = $tool->job(inputs => { ref => $reference, in => $input}, outputs => { out => $output});
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
