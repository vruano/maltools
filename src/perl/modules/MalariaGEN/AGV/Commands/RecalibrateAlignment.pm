package MalariaGEN::AGV::Commands::RecalibrateAlignment;

use strict;
use warnings;
use Moose;

use File::Copy qw(copy);

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);

use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long;
extends 'MalariaGEN::AGV::Command';

sub hidden {
  return 1;
}

sub help_summary {
   return 'realigns sample BAM files around indels';
}

sub help_text {
   return "recalibrate-alignment\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " recalibrate-alignment -i input.bam -o output.bam [-c count.bam] [ -r recalFile.csv ] [ -R fullReport.tgz ]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $input = undef;
  my $recalFile = undef;
  my $count_input;
  my $fullReport = undef;
  my $reference = undef;
  my $mask = reference_config->known_snps_mask_file();
  my $mask_format = reference_config->known_snps_mask_format();
  GetOptions("reference|ref=s" => \$reference, "i=s" => \$input, "o=s" => \$output, "c=s" => \$count_input, "recalFile|r=s" => \$recalFile, "fullReport|R=s" => \$fullReport);
  return $self->error_return("you must indicate an input bam to recalibrate") unless $input;
  return $self->error_return("you cannot use the standard input for the input bam") if ($input eq "-");
  -f $input or return $self->error_return("the input bam file '$input' does not exist or is not a regular file");
  return $self->error_return("you must indicate an output bam") unless $output;
  return $self->error_return("you cannot use the standard output for the output bam") if ($output eq "-");
  return $self->error_return("you cannot use the standard output for the recalfile") if ($recalFile && $recalFile eq "-");
  return $self->error_return("you cannot use the standard output for the full-report tarball") if ($recalFile && $fullReport eq "-");

  $count_input ||= $input;
  my $engine = $self->resolve_engine();
  my $ref = $reference || reference_config->file_name();
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-recalibrator');

  my $job = $tool->job(inputs => { ref => $ref, in => $input, ($count_input ne $input ? (recal_in => $count_input) : ()), mask => $mask, mask_format => $mask_format },
                       outputs => { out => $output, report => $fullReport, recal_file => $recalFile });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
  
}


1;
