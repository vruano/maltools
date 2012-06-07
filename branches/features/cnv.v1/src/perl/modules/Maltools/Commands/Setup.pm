package Maltools::Commands::Setup;

use Moose;

extends 'Maltools::Command';

use strict;
use warnings;

use Maltools::Config;
use Maltools::Tool;
use Getopt::Long;
use IO::File;

sub hidden { return 1; }

sub help_summary {
   return 'setups the environemnt for ' . $_[0]->cl_name . '';
}

sub help_text {
   my $prog = $_[0]->cl_name;
   my $summary = $_[0]->help_summary;
   return <<EOM
Command: 
     
 $prog - $summary

Synopsis:

 $prog setup 

Description:

 Perform all autoamatable set-up tasks so that $prog components can work propery.

EOM
;
}

sub execute {
  my ($self,$site,$params) = @_;
  print STDERR "checking c libraries...\n";
  return $self->error_return("missing libraries!") unless $self->setup_libraries;
  print STDERR "checking perl libraries...\n";
  return $self->error_return("missing perl modules!") unless $self->setup_perl;
  print STDERR "done!\n";
  return $self->ok_return;
}

sub setup_perl {
  my $self = shift;
#  my $config = Maltools::Config->new();
  my $config = $self->config();
  my $perl_deps_list_file = $config->get_resource('setup','perl-deps.list');
  defined $perl_deps_list_file or return $self->error_return('configuration does not define a perl dependency list');
  -f $perl_deps_list_file or return $self->error_return("could not find perl dependency list '$perl_deps_list_file'");
  open my $file , $perl_deps_list_file or return $self->error_return('could not open perl dependency list');
  my @missing_packages = ();
  while (my $package = <$file>) {
     chomp $package;
     next unless $package =~ /\S/;
     eval "require $package";
     $@ or next; # already there.
     push @missing_packages, $package;
  }
  my $package_list = join(" ",@missing_packages);
  return 1 unless @missing_packages;
  print STDERR <<EOT
Need to install the following missing perl packages " 

 $package_list.

Calling 'sudo cpan'...

EOT
;               
  system ("sudo cpan $package_list");
  $? and print STDERR "could not install or update perl packages required, sorry!\n", return 0;
  system ("Done installing perl dependencies!\n");
  return 1;
}


sub setup_libraries {
  my $self = shift;
  `locate libxml2.so`;
  return 1 unless $?;
  print STDERR "Missing libxml2, please install it\n";
  return 0;
}


1;
