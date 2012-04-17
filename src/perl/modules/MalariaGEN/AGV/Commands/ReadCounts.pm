package MalariaGEN::AGV::Commands::ReadCounts;

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
   return "read-counts\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " genotype -i input1.bam -i input2.bam ...  [-o output.vcf]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-"; 
  my @inputs = (); 
  my @regions  = ();
  my $threads = undef;
  GetOptions(
    "input|i=s@" => \@inputs,
    "output|o=s" => \$output,
    "region|r=s@" => \@regions,
    "threads|t=i" => \$threads,
  );
  $output ne "-" or return $self->error_return("you need to specify and output file name");
  $#inputs >= 0 or return $self->error_return("you need to speficy some samples");
  my $dc = data_config();
  my $rc = reference_config();
  my $ref = $rc->file_name;
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('gatk-readCounts');
  my $job = $tool->job(
      inputs => { 
          ref => $ref, 
          samples => \@inputs , 
          regions => [@regions],
          ref_uq => $rc->file_name(extension => '.uq'),
          ref_gff => $rc->file_name(extension => '.gff'),
          (defined $threads) ? ( threads => $threads) : () }, 
      outputs => { 
          out => $output  });

  my $engine = $self->resolve_engine;
  unless ($engine->run_job($job)) {
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotyping with message: " . $job->error_message);
  }
  else {
    return $self->ok_return;
  }
}



1;
