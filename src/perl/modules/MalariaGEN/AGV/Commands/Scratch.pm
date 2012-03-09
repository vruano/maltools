package MalariaGEN::AGV::Commands::Scratch;

use base 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(sanger_config);
use Getopt::Long;
use IO::File;
use base 'MalariaGEN::AGV::Command';

sub help_summary {
   return 'import-export files from data repository to Sanger\'s farm scratch area';
}

sub help_text {
   return "scratch:\n\t" .
          "import-exprot files from data repository to Sanger\'s farm scratch area\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " scratch [import|export] relative-file-or-dir\n" .
          "Subcommands:\n" .
          "\timport: copy from data repository to scratch\n";
          "\texport: copy from scratch to data repository\n"; 
}

sub execute {
  my ($self,$samtrak_site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $scommand = shift @ARGV;
  if ($scommand eq "import") {
    return $self->_import(@ARGV);
  }
  elsif ($scommand eq "export") {
    return $self->_export(@ARGV);
  }
  else {
    return $self->error_return("there is no such sub-command '$scommand'");
  }
}

sub _import {
  my ($self,@files) = @_;

  my $scratch = MalariaGEN::AGV::Engines::Scratch->new();
  foreach my $file (@files) {
    $scratch->import(file_or_dir => $file);
  }
  return $self->ok_return;
}

sub _export {
  my ($self,@files) = @_;

  my $scratch = MalariaGEN::AGV::Engines::Scratch->new();
  foreach my $file (@files) {
    $scratch->export(file_or_dir => $file);
  }
  return $self->ok_return;
}


1;
