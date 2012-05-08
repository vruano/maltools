package Maltools::Commands::CountCoverage;

use base 'Maltools::Command';

use strict;
use warnings;

use Maltools::Config qw(sanger_config reference_config jobs_config);
use Maltools::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'analyses the read depth coverage an a set of samples';
}

sub hidden {
   1
}

sub help_text {
   return "genotype\n\t" .
          $_[0]->help_summary . "\n" .
          "It generates a distribution description files in Json format containing the distribution all across samples and for each sample individually.\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " count-coverage -i input1.bam -i input2.bam -o output.json [-r REGION1 -r REGION2 ...] [-e errorfile.txt]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-"; my $errput = "-";
  my @inputs = ();
  my @regions = ();
  my $class = "realigned";
  my $fragments = 1;
  GetOptions("inputs|i=s@" => \@inputs, "fragments|f=i" => \$fragments, "o=s" => \$output,"e=s" => \$errput, "r=s@" => \@regions,"class|c=s" => \$class); 
  $output ne "-" or return $self->error_return("you need to specify and output file name");
  $#inputs >= 0 or return $self->error_return("you need to specify some input");
  my $rc = reference_config();
  my $ref = $rc->file_name;
  my $tool = Maltools::Tool->tool_by_name('gatk-countCoverage');
  my $job = $tool->job(
      inputs => { 
          ref => $ref, 
          samples => [@inputs] , 
          regions => [@regions],
          ref_gff => $rc->file_name(extension => '.gff'),
          ($fragments > 1) ? ( threads => $fragments) : () }, 
      outputs => { 
          out => $output });

  my $engine = Maltools::Engine->guess_engine;
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
