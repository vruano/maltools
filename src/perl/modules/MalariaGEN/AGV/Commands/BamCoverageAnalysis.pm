package MalariaGEN::AGV::Commands::BamCoverageAnalysis;

use strict;
use warnings;
use Moose;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long;
use File::Basename qw(fileparse);
use File::Spec::Functions qw(catfile);
use File::Path qw(remove_tree make_path);
use File::Copy qw(copy move);
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'calculate coverage statistics on bam files';
}

sub help_text {
   return "coverage-analysis\n\t" .
          $_[0]->help_summary . "\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $input = "-";
  my $region_size = 10000;
  my $remove_duplicates = 0;
  GetOptions("i=s" => \$input, "o=s" => \$output, "s=i" => \$region_size,"D" => \$remove_duplicates);
  $input eq "-" || -f $input or return $self->error_return("input file '$input' cannot be found or is not a regular file");  
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('bam-coverage');
  my $engine = MalariaGEN::AGV::Engines::Sanger->new();
  my $job = $tool->job(inputs => { in => $input, region_size => $region_size, remove_duplicates => $remove_duplicates },
                       outputs => { out => $output });
  if ($engine->run_job($job)) {
    return $self->ok_return;
  }
  else {
    return $self->error_return("Ban coverage analysis failed with exit code " . $job->return_code,$job->error_message);
  }
}


1;
