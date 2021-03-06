package Maltools::Vcf::Filters::loCVG;

use strict;
use warnings;
use Moose;

extends 'Maltools::Vcf::Filter';

has '+id' => ( default => 'loCVG' );
has '+description' => ( lazy => 1, default => sub { 'Low Coverage (thr=' . $_[0]->thr . ')' });
has 'thr' => (is => 'ro', default => '0.15');

sub passes {
 my ($self,$h) = @_; 
 my $info = $h->{INFO};
 my $ccp = $info->{CCP};
 return 1 unless defined $ccp;
 return $ccp >= $self->thr;
}



1;
