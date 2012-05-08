package Maltools::Commands::VariantInfo;

use strict;
use warnings;
use Moose;

use Maltools::Config qw(sanger_config reference_config jobs_config);
use Maltools::Engines::Sanger;
use Getopt::Long;
extends 'Maltools::Command';

sub help_summary {
   return 'add variant information to a vcf file';
}


sub hidden {
   return 1;
}

sub help_text {
   return "variant-info\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " variant-info INFO1 INFO2 INFO3 ... [-i input-file] [-o output-file]\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my $input = "-";
  GetOptions("i=s" => \$input, "o=s" => \$output);
   
  my $rc = reference_config();
  my $tool = Maltools::Tool->tool_by_name('vcf-info');
  my $engine = $self->resolve_engine();
  @ARGV = $self->translate_info_files($engine,@ARGV);  
  my $job = $tool->job(inputs => { in => $input, infos => [@ARGV], op => 'add' }, outputs => { out => $output});
  $engine->run_job($job);
  return $self->ok_return;
}


sub translate_info_files {
  my ($self,$engine,@ARGV) = @_;
  my @result = ();
  foreach my $info (@ARGV) {
    my %imports = ();
    while ($info =~ /=([^:]*)/g) {
      my $candidate = $1;
      if (exists($imports{$candidate})) {
        # nothing to do here.
      }
      elsif (-e $candidate) {
        $imports{$candidate} = $engine->import_resource($candidate);
      }
      else {
        # nothing to do here.
      }
    }
    foreach my $imported (keys %imports) {
        my $in_re = quotemeta($imported);
        my $out = $imports{$imported};
        $info =~ s/$in_re/$out/g;
    }
    push @result, $info;
  }
  print STDERR join(", ",@result),"\n";
  return @result;

}


1;
