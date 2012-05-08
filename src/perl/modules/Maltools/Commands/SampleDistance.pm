package Maltools::Commands::SampleDistance;

use strict;
use warnings;
use Moose;

use Maltools::Engines::Sanger;
use Maltools::Tool;
use Getopt::Long qw(:config no_ignore_case);
extends 'Maltools::Command';

sub help_summary {
   return 'create sample distance matrix ';
}

sub hidden {
  return 1;
}

sub help_text {
   return "sample-distance\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " sample-distance [-i INPUT-VCF1 -i INPUT-VCF2 ... | < INPUT-VCF]  [-o output-file | > output-file ]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my @inputs = ();
  my $input_file_list = undef;
  GetOptions("i|in|input=s" => \@inputs, "out|output|o=s" => \$output, "I=s" => \$input_file_list);
  if ($input_file_list) {
    my $fh = $input_file_list eq "-" ? \*STDIN : IO::File->new($input_file_list,'r') or die "cannot open input file list file '$input_file_list'";
    while (my $line = <$fh>) {
       chomp $line;
       next unless $line =~ /\S/;
       push @inputs,$line;
    }
    $fh->close();
  }
  @inputs = ("-") unless scalar(@inputs) > 0;
  my $tool = Maltools::Tool->tool_by_name('vcf-distanceMatrix');
  my $engine = $self->resolve_engine;
  my $job = $tool->job(inputs => { in => \@inputs },
                       outputs => { out => $output });
  my $result = $engine->run_job($job);
  return $self->ok_return if $result;
  return $self->error_return("error executing the vcf-distanceMatrix tool with return code " . $job->return_code,$job->error_message());
}


1;
