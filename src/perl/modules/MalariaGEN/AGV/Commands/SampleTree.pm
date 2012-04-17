package MalariaGEN::AGV::Commands::SampleTree;

use strict;
use warnings;
use Moose;

use MalariaGEN::AGV::Engines::Sanger;
use MalariaGEN::AGV::Tool;
use Getopt::Long qw(:config no_ignore_case);
use File::Copy qw(copy);
extends 'MalariaGEN::AGV::Command';

has '+engine_name' => ( default => 'local' );

sub help_summary {
   return 'create sample tree representations based on a set of distance matrices';
}
sub hidden {
  return 1;
}

sub help_text {
   return "sample-distance\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " sample-tree [-c|-C consensus-tree-file] [-T trees-file] [-i INPUT-DM1 -i INPUT-DM2 | < INPUT-DM] [-o output-file | > output-file ]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my @inputs = ();
  my $output_consensus = 0;
  my $consensus_file = undef;
  my $trees_file = undef;
  GetOptions("i|in|input=s" => \@inputs, "out|output|o=s" => \$output, "c!" => \$output_consensus, "T|trees-file=s" => \$trees_file, "C|consensus-file=s" => \$consensus_file);
  
  print STDERR "Warning: does not make much sense to combine -c -T and -C all at the same time\n" if ($output_consensus && $trees_file && $consensus_file); 
  $trees_file ||= $output_consensus ? undef : ($output || "-");
  $consensus_file ||= $output_consensus ? ($output || "-") : undef;
  @inputs = ("-") unless scalar(@inputs) > 0;
  my @temporal_files = ();
  my $engine = $self->resolve_engine;

  # Merge all input distance matrices into a single file
  my $input = $self->merge_inputs($engine,\@temporal_files,@inputs);

  # Decide where to store the bootstrapped trees:
  my $trees_temp_file;
  if ($trees_file && $trees_file ne "-") {
    $trees_temp_file = $trees_file;
  }
  else {
    $trees_temp_file = $engine->tempfile();
    push @temporal_files, $trees_temp_file;
  }

  # Run neighbor:
  
  my $neighbor = MalariaGEN::AGV::Tool->tool_by_name("phylip-neighbor");
  my $neighbor_job = $neighbor->job(inputs => { in => $input,
                                                  format  => 'lower-triangular',
                                                  seed => 13 },
                                      outputs => { out_tree => $trees_temp_file, out_file => undef });
  my $result = $engine->run_job($neighbor_job);
  unless ($result) {
    foreach my $f (@temporal_files) { unlink $f };
    return $self->error_return("error executing the phylip-neighbor tool with return code " . $neighbor_job->return_code,$neighbor_job->error_message());
  }
  if ($trees_file && $trees_file ne $trees_temp_file) {
    my $fh = $trees_file eq "-" ? \*STDOUT : IO::File->new($trees_file,"w");
    copy($trees_temp_file,$fh);
    $fh->close;  
  } 

  # Finish here if the consensus tree was not requested.
  unless ($consensus_file) {
    foreach my $f (@temporal_files) { unlink $f };
    return $self->ok_return();   
  } 
   
  # Run consense:
  my $consense = MalariaGEN::AGV::Tool->tool_by_name("phylip-consense");
  my $consense_job = $consense->job(inputs => { in => $trees_temp_file } , outputs => { out_tree => $consensus_file, out_file => undef });
  $result = $engine->run_job($consense_job);
  foreach my $f (@temporal_files) { unlink $f };
  unless ($result) {
    return $self->error_return("error executing the phylip-consense tool with return code " . $consense_job->return_code,$consense_job->error_message());
  }
  return $self->ok_return;
}

sub merge_inputs {
  my ($self,$engine,$temporal_files,@inputs) = @_;
  if (scalar(@inputs) > 1) {
    my $convined_input = $engine->tempfile();
    push @$temporal_files, $convined_input;
    my $cifh = IO::File->new($convined_input,'w')  
         or die "could not open combined input file '$convined_input' for writing";
    foreach my $input (@inputs) {
      my $ifh = $input eq "-" ? \*STDIN : IO::File->new($input,'r');
      copy($ifh,$cifh) 
         or die "could not copy input content in '$input' into the combined file '$convined_input'";
      $ifh->close() unless $input eq "-";
    }
    $cifh->close();
    return $convined_input;
  }
  else {
    return $inputs[0];
  }
}


1;
