package Maltools::Commands::CrossCalls;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);

use Maltools::Config qw(data_config sanger_config reference_config jobs_config);

use Maltools::Engines::Sanger;
use Getopt::Long qw(:config no_ignore_case);
extends 'Maltools::Command';

sub hidden {
  return 1;
}

sub help_summary {
   return 'Add cross and genotype and error calls to existing vcf';
}



sub help_text {
   return "cross-calls\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " cross-call -i input.vcf -p Parent1 -p Parent2 -o output.vcf --matrices matrices-dir --report report-dir\n"; 
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $matrices = undef;
  my $report = undef;
  my $input = undef;
  my @parents = ();
  GetOptions("matrices|m=s" => \$matrices, "report|r=s" => \$report, 
      "parent|p=s@" => \@parents, "output|out|o=s" => \$output,  "input|i=s" => \$input );

  $#parents == 1 or return $self->error_return("you must indicate exactly two parents");
  $input && -f $input or return $self->error_return("the input file '$input' cannot be reached or is not a regular file");
  $matrices && -d $matrices or return $self->error_return("the before file '$matrices' does not exists");
  $report && -d $report or return $self->error_return("the report directory '$report' does not exists or is not a regular file");
  $output && ( !-e $output || -d $output) or return $self->error_return("the output '$output' must be provided, not exist or if exists be a directory");
   
  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('crossCalls');
  my 
  $job = $tool->job(inputs => { in => $input, parents =>\@parents, matrices => $matrices, report => $report  }, outputs => { out => $output});
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
