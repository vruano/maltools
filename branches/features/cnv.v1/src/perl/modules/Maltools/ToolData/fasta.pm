package Maltools::ToolData::fasta;

use strict;
use warnings;
use Moose;

extends 'Maltools::ToolData::file';

has 'indexed' => ( is => 'ro', default => 0);
has '+has_files' => ( default => 1 );

sub file_names {
  my ($self,$value) = @_;
  my @result = $self->SUPER::file_names($value);
  if ($self->indexed) {
    push @result, (ref($value) eq "ARRAY" ? @$value : ( $value . ".fai" ));
  }
  return wantarray ? @result : \@result;
}

1;
