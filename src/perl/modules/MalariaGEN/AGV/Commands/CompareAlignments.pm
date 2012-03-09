package MalariaGEN::AGV::Commands::CompareAlignments;

use base 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'compare two alignments';
}

sub help_text {
   return "genotype\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " compare-alignments -l left.bam -r right.bam -o output.diff\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $left = undef;
  my $right = undef;
  GetOptions(
    "left|l=s" => \$left,
    "rigth|r=s" => \$right,
    "ouput|o=s" => \$output,
  );
 
  $left && -f $left or return $self->error_return("the left alignment file cannot be accessed or is not a regular file");
  $right && -f $right or return $self->error_return("the right alignment file cannot be accessed or is not a regular file");
  $output or return $self->error_return("you need to specify a output file");
  
  my $rc = reference_config();
  my $ref = $rc->file_name;
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-compareAlignments'); 
  my $job = $tool->job(
      inputs => { 
          ref => $ref, 
          left => $left,
          right => $right,
      },
      outputs => { 
          out => $output 
  });

  my $engine = MalariaGEN::AGV::Engine->guess_engine;
  unless ($engine->run_job($job)) {
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotyping with message: " . $job->error_message);
  }
  else {
    return $self->ok_return;
  }
}

1;
