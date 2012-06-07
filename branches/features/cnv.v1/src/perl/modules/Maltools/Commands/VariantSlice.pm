package Maltools::Commands::VariantSlice;

use strict;
use warnings;
use Moose;

use Maltools::Config qw(data_config 
       sanger_config reference_config jobs_config);
use Maltools::Engines::Sanger;
use Getopt::Long qw(:config pass_through no_getopt_compat no_ignore_case);
extends 'Maltools::Command';

sub help_summary {
   return 'slice out variants based on whether they passed filters or not';
}


sub hidden {
   return 1;
}

sub help_text {
   return "variant-slice\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " variant-slice [-m (first|last)] [-S include_sample_file|-s SAMPLE1 -s SAMPLE2 ...] [condition1 condition2 ...]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my $input = "-";
  my $mode = "last";
  my $default = 'exclude';
  my @include_filters = ();
  my @exclude_filters = ();
  my @samples = ();
  my $sample_file = undef;
  GetOptions("i=s" => \$input, "o=s" => \$output, "m=s" => \$mode, "d=s" => \$default, "s=s@" => \@samples, "S=s" => \$sample_file);
  $mode = lc($mode); $default = lc($default);
  die "illegal logical mode '$mode' must be either 'first' or 'last'" unless $mode eq 'first' || $mode eq 'last';
  die "illegal default action '$default' must be either 'include' or 'exclude'" unless grep { $_ eq $default } qw(include exclude inc exc incl excl - +);
  $default = index($default,"in") == 0 || $default eq "+" ? "include" : "exclude";
  my %samples = map { $_ => 1 } @samples;
  if ($sample_file) {
    my $fh = IO::File->new($sample_file) or die "could not open sample list file '$sample_file' for reading";
    while (my $line = <$fh>) {
      chomp $line;
      my @s = split(/\s+/,$line);
      foreach my $s (@s) { $samples{$s} = 1 }; 
    }
    $fh->close();
  }
  @samples = sort keys %samples;
  my $tool = Maltools::Tool->tool_by_name('vcf-slice');
  my $engine = $self->resolve_engine();
  my $job = $tool->job(inputs => { in => $input, filters => \@ARGV, mode => $mode, default => $default, samples => [@samples] }, outputs => { out => $output});
  my $result = $engine->run_job($job);
  return $self->ok_return if $result;
  return $self->error_return("error executing tool vcf-slice with return code " . $job->return_code,$job->error_message());
}


1;
