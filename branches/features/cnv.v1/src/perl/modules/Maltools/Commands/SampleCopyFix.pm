package Maltools::Commands::SampleCopyFix;

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
   return 'copy and alignment fixing existing issues with @PG headers';
}

sub help_text {
   return "sample-copy-fix\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " sample-copy-fix -i input-bam -o output-file\n"; 
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $input = undef;
  GetOptions("input|in|i=s" => \$input, "output|out|o=s" => \$output);
   
  -f $input or return $self->error_return("the input alignment '$input' is unreachable or is not a regular file");
  $output or return $self->error_return("the output alignment '$output' must be specified");
   
  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('sampleCopyFix');
  my $job = $tool->job(inputs => { in => $input }, outputs => { out=> $output });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
