package MalariaGEN::AGV::ROD::UniquenessScore;

use strict;
use warnings;

use Moose;

extends 'MalariaGEN::AGV::ROD';

sub _process_line {
  my $self = shift;
  my $line = shift;
  my $prev = shift;
  return () if $line =~ /^\s*#/;
  return () if $line !~ /\S/;
  my ($seq,$pos,$score) = split(/\t/,$line);
  return ($seq,$pos,Feature->new( seq => $seq, pos => $pos, score => $score));
}

1;
