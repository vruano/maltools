package Maltools::Commands::CoverageAnalysis;

use strict;
use warnings;
use Moose;

use Maltools::Config qw(data_config 
       sanger_config reference_config jobs_config);
use Maltools::Engines::Sanger;
use Getopt::Long;
use File::Basename qw(fileparse);
use File::Spec::Functions qw(catfile);
use File::Path qw(remove_tree make_path);
use File::Copy qw(copy move);
extends 'Maltools::Command';

sub help_summary {
   return 'create counter group files to reflect coverage depth in variants';
}

sub hidden {
  return 1;
}

sub help_text {
   return "coverage-analysis\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " coverage-analysis [-g GFF annotation-file] [-i input-file] [-o output-directory] \n\n" .
          "if not output directory is provided output files names will start by input-file's basename if\n" .
          "provided if not the empty prefix '' will be used instead.\n" -
          "Output file name patterns follow:\n\n" .
          "\t<prefix>-cvg-<seq>.json for each sequences or chromosome includes the counter group structure in JSON format\n" .
          "\t<prefix>-cvg-all.json       the overall counter group\n" .
          "\t<prefix>-cvg-nuclear.json   only those nuclear sequences (i.e. non-organelle DNA included)\n" .
          "\t<prefix>-cvg-[<seq>|all|nuclear].png the coverage density plots\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = undef;
  my $input = "-";
  my $gff = undef;
  my $skip_filtered = 0;
  GetOptions("i=s" => \$input, "o=s" => \$output, "g=s" => \$gff,"F+" => \$skip_filtered);
  $skip_filtered = 1 if $skip_filtered;
  
  unless (defined $gff) {
   $gff = reference_config()->file_name(extension => '.gff');
  }  
  -f $gff or die "gff annotation file $gff does not exists";
 
  my $tool = Maltools::Tool->tool_by_name('vcf-coverage');
  my $engine = Maltools::Engines::Sanger->new();
  my $out_dir = $output || $engine->tempdir();
  -d $out_dir || make_path($out_dir) or die "$out_dir is not a directory";
  my $job = $tool->job(inputs => { gff => $gff, in => $input, skip_filtered => $skip_filtered }, outputs => { out => $out_dir });
  $engine->run_job($job);
  # copy the output files over
  unless ($output) {
    opendir(OUTDIR,$out_dir);
    my $prefix = $input eq "-" ? "" : fileparse($input,qr/\.[^.]*/) . "-";
    while(my $file = readdir(OUTDIR)) { 
       next if $file =~ /^\./;
       #next unless -f $file;
       print STDERR "the move from " . catfile($out_dir,$file) . " to " . $prefix . $file . "\n";
       move(catfile($out_dir,$file),$prefix . $file);
    } 
    closedir(OUTDIR);
    remove_tree($out_dir);
  }
  return $self->ok_return;
}


1;
