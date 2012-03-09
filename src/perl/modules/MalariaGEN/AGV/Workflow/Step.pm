package MalariaGEN::AGV::WorkflowStep;

use strict;
use warnings;

use Moose;

our $anonymous_step_id = "ws_00000";

has 'tool' => ( is => 'ro', isa => 'Str', required => 1);
has 'inputs' => ( is => 'ro', isa => 'HashRef', default => sub {{}});
has 'outputs' => ( is => 'ro', isa => 'HashRef', default => sub {{}});
has 'id' => ( is => 'ro', isa => 'Str', default => sub { return $anonymous_step_id++; });
has 'num' => ( is => 'ro', writer => '_set_num', default => sub { -1 } );


around BUILDARGS => sub {
  my $orig = shift;
  my $class = shift;
  my %args = @_;
  $args{tool} = MalariaGEN::AGV::Tool->tool_by_name($args{tool}) if defined $args{tool} && ref($args{tool}) eq "";
  return $class->$orig(%args);
};


sub names {
  my $self = shift;
  return wantarray ? ($self->id, $self->num) : [$self->id, $self->num];
}

1;
