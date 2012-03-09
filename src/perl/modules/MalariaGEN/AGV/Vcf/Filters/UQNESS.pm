package MalariaGEN::AGV::Vcf::Filters::UQNESS;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::Vcf::Filter';

has '+id' => ( default => 'UQNESS' );
has '+description' => ( default => 'Uniqueness score threshold filter');
has 'thr' => (is => 'ro', default => '26');
has 'ethr' =>  ( is => 'ro', lazy => 1, builder => '_ethr_builder' );

sub _ethr_builder {
 my ($self);
 my $thr => $self->thr;
 return $thr if $thr <= 2;
 return 2 + (log($thr - 2));
}

sub passes {
 my ($self,$h) = @_; 
 my $info = $h->{INFO};
 my $uq = $info->{UQ};
 return 1 unless defined $uq;
 my $thr = $self->ethr;
 return 1 if $uq <= $thr;
 my $gtypes = $h->{gtypes};
 foreach my $sample (keys %$gtypes) {
   my $gt = $gtypes->{$sample}->{GT};
   next unless defined $gt && index($gt,'.') == -1;
   my $a1 = substr($gt,0,1);
   my $a2 = substr($gt,2,1);
   return 0 unless $a1 eq $a2;
 }
 return $uq;
}

1;
