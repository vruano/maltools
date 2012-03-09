package MalariaGEN::AGV::Commands::RecalAnalysis;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);

use MalariaGEN::AGV::Config qw(data_config sanger_config reference_config jobs_config);

use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long qw(:config no_ignore_case);
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'perform the recalibration quality score analysis';
}

sub help_text {
   return "recal-analysis\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " recal-analysis -s sample -v variant-calls.vcf -b before.bam -a after.bam -o output-directory [--reference reference]?\n"; 
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $before = undef;
  my $after = undef;
  my $vcalls = undef;
  my $sample = undef;
  my $baq = 1;
  my $reference = undef;
  GetOptions("before|b=s" => \$before, "after|a=s" => \$after, 
      "sample|s=s" => \$sample, "output|out|o=s" => \$output, 
      "reference|ref|r=s" => \$reference, "vcf|v=s" => \$vcalls,
             "BAQ|B!" => \$baq );
  $vcalls && -f $vcalls or return $self->error_return("the variant calls file '$vcalls' does not exists");
  $before && -f $before or return $self->error_return("the before file '$before' does not exists");
  $after && -f $after or return $self->error_return("the after fiel '$after' does not exists or is not a regular file");
  $sample or return $self->error_return("you need to indicate an input sample");
  $output && ( !-e $output || -d $output) or return $self->error_return("the output '$output' must be provided, not exist or if exists be a directory");
  my $rc = reference_config();
   $reference ||= $rc->file_name();
  -f $reference or return $self->error_return("there reference '$reference' could not be reached");
   
  my $engine = $self->resolve_engine();
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('recalAnalysis');
  my 
  $job = $tool->job(inputs => { ref => $reference, before => $before, after => $after, 
          sample => $sample, baq => $baq, vcf => $vcalls, }, outputs => { out_dir => $output});
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
