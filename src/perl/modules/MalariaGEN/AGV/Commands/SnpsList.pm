package MalariaGEN::AGV::Commands::SnpsList;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);

use MalariaGEN::AGV::Config qw(data_config 
  sanger_config reference_config jobs_config);

use MalariaGEN::AGV::Engines::Sanger;
use Getopt::Long;
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'filters an alignment reads using the Snp-o-Matic approach (Bwa2Som)';
}

sub help_text {
   return "snps-list\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " snps-list -i input-bam -o output.snps [--reference reference]? [--snpList alternativeCandidateList]?\n"; 
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $input = undef;
  my $reference = undef;
  my $snpList = undef;
  my $minq = undef;
  my $rc = reference_config();
  GetOptions("input|in|i=s" => \$input, "output|out|o=s" => \$output, "reference|r=s" => \$reference,
             "snpList|snps=s" => \$snpList, "minq=i" => \$minq);
  $snpList ||= $rc->possible_snps_list();
  print STDERR "$snpList\n";
  $reference ||= $rc->file_name();
  $minq = 27 unless defined $minq;
   
 
  -f $input or return $self->error_return("the input alignment '$input' is unreachable or is not a regular file");
  -f $snpList or return $self->error_return("the possible snp-list '$snpList' is unreachable or is not a regular file");
  -f $reference or return $self->error_return("the reference '$reference' is unreachable or is not a regular file");
   
  my $engine = $self->resolve_engine();
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('bam2Snps');
  my $job = $tool->job(inputs => { ref => $reference, in => $input, candidates => $snpList, minq => $minq}, outputs => { out => $output});
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during realignment tool execution",$job->error_message());
  }
  return $self->ok_return;
}


1;
