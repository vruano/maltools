package MalariaGEN::AGV::Commands::VariantBootstrap;

use strict;
use warnings;
use Moose;

use MalariaGEN::AGV::Engines::Sanger;
use MalariaGEN::AGV::Tool;
use Getopt::Long qw(:config no_ignore_case);
extends 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'bootstrap vcf file';
}

sub help_text {
   return "variant-bootstrap\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " variant-bootstrap [-i input-vcf| < input-vcf] -o output-dir [-S summary_file] [-c chromose] -r num-repeats -s size-per-repeat [--seed random_seed] [-p file-pattern] \n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $out_dir = undef;
  my $input = "-";
  my $repeats = 100;
  my $size = 10000;
  my $seed = undef;
  my $pattern = "bootstrap_%04d.vcf";
  my $summary_file = undef;
  my $chromosome = undef;
  GetOptions("chromosome|chr|c=s" => \$chromosome, "in|input|i=s" => \$input, "out|out-dir|o=s" => \$out_dir,
             "size|s=i" => \$size, "repeats|r=i" => \$repeats, "pattern|p=s" => \$pattern, "seed|S=i" => \$seed,
             "summary-file|S=s" => \$summary_file);
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('vcf-bootstrap');
  my $engine = $self->resolve_engine;
  my $job = $tool->job(inputs => { in => $input, 
                                   seed => $seed, 
                                   repeats => $repeats,
                                   size => $size,
                                   pattern => $pattern,
                                   (defined $chromosome ? (chr => $chromosome) : ()),
                                 },
                       outputs => { out_dir => $out_dir,
                                    summary => $summary_file });
  my $result = $engine->run_job($job);
  return $self->ok_return if $result;
  return $self->error_return("error executing tool vcf-slice with return code " . $job->return_code,$job->error_message());
}


1;
