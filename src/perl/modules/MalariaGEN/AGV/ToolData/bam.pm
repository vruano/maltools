package MalariaGEN::AGV::ToolData::bam;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::ToolData::file';

has 'indexed' => ( is => 'ro', default => sub { 0 } );

sub file_names {
  my ($self,$value) = @_;
  my @result = ref($value) eq "ARRAY" ?  @$value : ( $value );
  if ($self->indexed) {
    push @result, ((ref($value) eq "ARRAY") ? (map { $_ . ".bai" } @$value) : ( $value . ".bai" ));
  }
  return wantarray ? @result : \@result;
}

1;
