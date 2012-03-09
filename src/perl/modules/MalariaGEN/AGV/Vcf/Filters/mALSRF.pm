package MalariaGEN::AGV::Vcf::Filters::mALSRF;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::Vcf::Filter';

has '+id' => ( default => 'mALSRF' );
has '+description' => ( lazy => 1, 
    default => sub { 'Minor Allele Sample Read Frequency filter. Filters out positions where any of the two alleles (but specially the minor allele) is seen in any sample in more than a minimum threshold number of reads (thr=' . $_[0]->thr . ')' } );
has 'thr' => (is =>  'ro' , default => 10);

sub passes {
 my ($self,$h) = @_; 

 my $alt = $h->{ALT};
 my $gt = $h->{gtypes};
 # We do not check on triallelic or not variant positions.
 my @allele_maxs = map { 0 } (0 .. scalar(@$alt));
 foreach my $sample (keys %$gt) {
   my $counts = $gt->{$sample}->{AD} or next;
   my @counts = split(/,/,$counts);
   for (my $i = 0; $i < scalar(@counts); $i++) { 
      $allele_maxs[$i] = $counts[$i] if $allele_maxs[$i] < $counts[$i];
   } 
 }
 # Remove totally missing variants.
 @allele_maxs = grep { $_ } @allele_maxs;
 # Pass non-biallelic variants, we do not care about them.
 return 1 unless scalar(@allele_maxs) != 2;
 my $thr = $self->thr;
 return scalar(grep { $_ >= $thr } @allele_maxs) == 2;
}



1;
