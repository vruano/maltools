package Maltools::ToolData::reference;

use strict;
use warnings;
use Moose;

extends 'Maltools::ToolData::fasta';

has 'indexed' => ( is => 'ro', default => 1);
has 'dicted' => (is => 'ro', default => 1);
has 'sam_hdr' => (is => 'ro', default => 0);
has '+has_files' => ( default => 1 );

sub file_names {
  my ($self,$value) = @_;
  my @result = $self->SUPER::file_names($value);
  if ($self->indexed) {
    push @result, (ref($value) eq "ARRAY" ? (map { $_ . ".fai" } @$value) : ( $value . ".fai" ));
  }
  if ($self->dicted) {
    push @result, (ref($value) eq "ARRAY" ? (map { $_ . ".dict" } @$value) : ($value . ".dict"));
  }
  if ($self->sam_hdr) {
    push @result, (ref($value) eq "ARRAY" ? (map { $_ . ".sam-hdr" } @$value) : ($value . ".sam-hdr"));
  }
  return wantarray ? @result : \@result;
}

1;
