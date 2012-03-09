package MalariaGEN::AGV::Workflow::Var;

use strict;
use warnings;

use Moose;

has 'name' => (is => 'ro', isa => 'Str', default => sub { return $_->names->[0] if exists($_->{names}) });
has 'step' => (is => 'ro');
has 'type' => (is => 'ro', isa => 'Str', required => 1);

sub names {
  my $step = $_->step;
  my $name = $_->name;
  return wantarray ? ($name) : [$name] unless defined $step;
  my @result = map { $name . "::" . $_ } $step->names;
  return wantarray ? @result : \@result;
}


1;
