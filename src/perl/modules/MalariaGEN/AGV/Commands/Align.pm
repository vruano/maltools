package MalariaGEN::AGV::Commands::Align;

use base 'MalariaGEN::AGV::Command';

use Cwd qw(realpath);
use strict;
use warnings;
use MalariaGEN::AGV::Config;
use MalariaGEN::AGV::Tool;
use MalariaGEN::AGV::Alignment;
use Getopt::Long;
use IO::File;
use File::Basename qw(dirname);

sub help_summary {
   return 'align a number of sequence files in bam format against a  or the reference';
}

sub help_text {
   my $summary = $_[0]->help_summary;
   my $cl_name = $_[0]->cl_name;
   return <<EOM;
   
Command:
   
   align - $summary.

Syntaxis:
      
   $cl_name align -i input.bam1 -o output.bam [-r reference.fasta ] [--force]
         
      aligns the content of the input files in a new output bam file.
      If the input files are already aligned against the reference provided it only will
      merge their content unless you use --force.
   
      Thinks are required to this program to consider a file to be aligned against the
      referernce:
  
      1. the input header must contain a reference dictionary matching the one
         of the target reference: name of chromosomes, order and md5 check sums must be the same.
      2. The alignment algorithm in the header must be the same as the targeted one.

      If the input files are all aligned against the reference and you do not force a realignemt
      this will be in effect be the same as a bam merge.

EOM
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my @inputs = ();
  my $output = "-";
  my $force = 0;
  my $reference = 0;
  GetOptions("i|input=s@" => \@inputs,"o|output=s" => \$output, "r|reference=s" => \$reference,
             "f|force!" => \$force);

  $reference ||= MalariaGEN::AGV::Config->new()->get_reference_sequences();

  return $self->error_return("you must indicate  at least one input ") unless $#inputs >= 0;
  my @bad_inputs = grep { ! -e $_ } @inputs;
  return $self->error_return("some input files are missing or cannot be reached: " . join(" ", @bad_inputs)) if $#bad_inputs >= 0;
  @bad_inputs = grep { ! -f $_ } @inputs;
  return $self->error_return("some input files are not regular files: " . join(" ",@bad_inputs)) if $#bad_inputs >= 0;
  @bad_inputs = grep { ! -r $_ } @inputs;
  return $self->error_return("some input files are not regular files: " . join(" ",@bad_inputs)) if $#bad_inputs >= 0;
  
  return $self->error_return("you must indicate and output file") unless $output;
  return $self->error_return("the output file exists and is not a regular file") if -e $output && !-f $output;

  my $out_dir = dirname($output);
  return $self->error_return("the output directory '$out_dir' does not exists or is not a directory") unless -e $out_dir && -d $out_dir;
  return $self->error_return("you cannot write in the output directory '$out_dir'") unless -w $out_dir;

  my $tool;
  my $job;
  unless ($force) {
    my $in_aln = MalariaGEN::AGV::Alignment->new(file => $inputs[0]);
    my $in_ref = MalariaGEN::AGV::Reference->new(sequence_file => $reference);
    if ($in_aln->valid_reference($in_ref)) {
       if ($in_aln->order ne 'coordinate') {
          $tool = MalariaGEN::AGV::Tool->tool_by_name("picard-sortSam");
          $job = $tool->job(inputs => { in => $inputs[0] }, outputs => { out => $output });
       }
       else {
          return $self->ok_return;
       }
    }
  }
  unless ($tool) {
    $tool = MalariaGEN::AGV::Tool->tool_by_name("align");
    $job = $tool->job(inputs => { in => $inputs[0], force => $force, ref => $reference },
                        outputs => { out => $output });
  }
  my $engine = $self->resolve_engine;
  $engine->run_job($job) 
	or $self->error_return("error executing the align command with exit code : " . $job->return_code,$job->error_message );
  return $self->ok_return;
}

1;
