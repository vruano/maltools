package Maltools::ToolData::file;

use strict;
use warnings;
use Moose;

extends 'Maltools::ToolData::Type';

has '+has_files' => ( default => 1 );

sub file_names {
  my ($self,$value) = @_;
  return wantarray ? @$value : $value if ref($value) eq "ARRAY";
  return wantarray ? ($value) : [$value];
}

1;
