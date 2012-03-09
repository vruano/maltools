package MalariaGEN::AGV::Vcf::Filters::hiCVG;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::Vcf::Filter';

has '+id' => ( default => 'hiCVG' );
has '+description' => ( lazy => 1, default => sub { 'High Coverage (thr=' . $_[0]->thr . ')' });
has 'thr' => (is => 'ro', default => '0.85');

sub passes {
 my ($self,$h) = @_; 
 my $info = $h->{INFO};
 my $ccp = $info->{CCP};
 return 1 unless defined $ccp;
 return $ccp <= $self->thr;
}



1;
