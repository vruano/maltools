package MalariaGEN::AGV::Commands::BamSetAnalysis;

use base 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(data_config 
       sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub hidden {
  return 1;
}

sub help_summary {
   return 'genotypes a set of samples';
}

sub help_text {
   return "bam-set-analysis\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " -i input.bam1 -i input.bam2 ... [--cvg-out\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my @regions = ();
  my $fragments = 1;
  my $bq_output = undef;
  my $baq_bq_output = undef;
  my $rbsq_output = undef;
  my $rmq_output = undef;
  my $cvg_output = undef;
  my $output_base = undef;
  my @files = ();
  GetOptions(
    "output|out|o=s" => \$output_base,
    "input|i=s@" => \@files,
    "fragments|f=i" => \$fragments, 
    "r=s@" => \@regions,
    "bq-output|bq=s" => \$bq_output,
    "rbsq-output|rq=s" => \$rbsq_output,
    "rmq-output|rm=s" => \$rmq_output,
    "baq-bq-output|baqBq=s" => \$baq_bq_output,
    "cvg-output|cvg=s" => \$cvg_output); 
  
  $#files >= 0 or return $self->error_return("you must indicate at least one input alignment");
  grep { $_ } ($output_base,$cvg_output,$rmq_output,$rbsq_output,$baq_bq_output,$bq_output) or return $self->error_return("you must indicate some analysis output");
  if ($output_base) {
    $cvg_output ||= $output_base . '.cvg.json';
    $bq_output ||= $output_base . '.baseq.json';
    $rbsq_output ||= $output_base . '.rbsq.json';
    $rmq_output ||= $output_base . '.mapq.json';
    $baq_bq_output ||= $output_base . '.baqbq.json';
  }
  my $dc = data_config();
  my $rc = reference_config();
  grep { ! -f $_ } @files and return $self->error_return("some of the input bam provided does not seem to exist!!!: " . join(", ", grep { ! -f $_ } @files));
  my $ref = $rc->file_name;
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-statistics');
  my $job = $tool->job(
      inputs => { 
          ref => $ref, 
          samples => [@files] , 
          regions => [@regions],
          interval_list => $rc->file_name(extension => '.interval_list'),
          ref_uq => $rc->file_name(extension => '.uq'),
          ref_gff => $rc->file_name(extension => '.gff'),
          ($fragments > 1) ? ( threads => $fragments) : () }, 
      outputs => { 
          $bq_output ? (bq_out => $bq_output):(),
          $baq_bq_output ? (baq_bq_out => $baq_bq_output):(),
          $rbsq_output ? (rbsq_out => $rbsq_output):(),
          $rmq_output ? (rmq_out => $rmq_output):(),
          $cvg_output ? (cvg_out => $cvg_output):() }
  );

  my $engine = $self->resolve_engine;
  unless ($engine->run_job($job)) {
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotyping with message: " . $job->error_message);
  }
  else {
    return $self->ok_return;
  }
}

sub sample_to_file {
  my ($self,$file,$class,$dc) = @_;
  if ($class eq "file") {
    return undef unless -f $file;
    return undef unless -f $file . ".bai";
  }
  return $dc->sample_bam(ox_code => $file, indexed => 1, class => $class);
}


1;
