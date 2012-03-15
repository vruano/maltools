package MalariaGEN::AGV::Command;

use strict;
use warnings;
use Scalar::Util qw(blessed);
use Moose;
use MalariaGEN::AGV::Engine;
use MalariaGEN::AGV::Manifest::Config;
use MalariaGEN::AGV::Config;
use File::Basename qw(basename);

has 'engine_name' => ( is => 'rw', isa => 'Str',  lazy => 1, default => 'auto' );
has 'procs' => (is => 'rw', isa => 'Num', lazy => 1, default => 1 );
has 'running_options' => (is => 'rw', isa => 'HashRef', default => sub { {} });
has 'samtrak' => (is => 'rw');
has 'config' => (is => 'ro', isa => 'MalariaGEN::AGV::Manifest::Config', lazy => 1, builder => '_build_config' );


sub resolve_engine {
  my ($self) = @_;
  my $result = MalariaGEN::AGV::Engine->engine_by_name($self->engine_name);
  $result->running_options($self->running_options);
  return $result;
}

sub hidden {
  return 0;
}

sub cl_name {
  return basename($0);
}

sub command_class_by_name {
  my ($class,$command_name) = @_;
  my @cnparts = split(/\-+/,$command_name);
  $command_name = 'MalariaGEN::AGV::Commands::' . join("",map { uc(substr($_,0,1)) . substr($_,1) } @cnparts);
  eval "require $command_name";
  $@ and return undef;
  return $command_name;
}

sub instance {
  my $cls = $_[0]->command_class_by_name($_[1]) or die "there is not such a command $_[1]\n";
  return $cls->new();
}

sub command_by_name {
  my $cls = $_[0]->command_class_by_name($_[1]) or die "there is not such a command $_[1]\n";
  return $cls->new();
}

sub help {
  my ($self,$format) = @_;
  my $help_method = $self->can('help_' . $format) || $self->can('help_default');
  return $help_method->($self,$format);
}

sub help_summary {
  return '';
}

sub help_default {
  my ($self,$format) = @_;
  return 'unknown help format ' . $format;
}

sub help_text {
  my ($self) = @_;
  return 'no help available';
}

sub help_wiki {
  my ($self) = @_;
  my $result = "==== " . $self->name . " ====\n";
  $result .= $self->help_text . "\n";
  return $result; 
}

sub name {
  my ($self) = @_;
  my $class = blessed($self) ? ref($self) : $self;
  $class =~ /MalariaGEN::AGV::Commands::(.+)$/;
  my @parts = split(/(?<=[a-z0-9])(?=[A-Z])/, $1);
  @parts = map { $_ =~ /^[A-Z][A-Z]*$/ ? $_ : lc($_) } @parts;
  my $result = join("-",@parts);
  return $result;
}

sub run {
  my ($self,$site,@args) = @_;
  $self->samtrak($site);
  return $self->execute($site,@args);
}

sub execute {
  my ($self,$site,@args) = @_;
  return $self->ok_return();
}

sub ok_return {
  return 0;
}

sub error_return {
  my ($self,$error,$details) = @_;
  print STDERR "Error occurred with message: $error\n";
  print STDERR "Details follow:\n\t" . $details . "\n" if defined $details;
  return 1;
}

sub _build_config {
  my $self = shift;
  return MalariaGEN::AGV::Config->new();
}

1;
