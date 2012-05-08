package Maltools::Vcf::Filters::mALF;

use strict;
use warnings;
use Moose;

extends 'Maltools::Vcf::Filter';

has '+id' => ( default => 'mALF' );
has '+description' => ( lazy => 1, default => sub { 'Minor Allele Low Frequency filter. Filters out positions where the minor allele (that with the lower number of total reads) does not appear in at least a threshold fraction of reads (thr=' . $_[0]->thr . '%)' } );
has 'thr' => (is =>  'ro' , default => 0.01);

sub passes {
 my ($self,$h) = @_; 

 my $alt = $h->{ALT};
 my $gt = $h->{gtypes};
 # We do not check on triallelic or not variant positions.
 my $total = 0;
 my @allele_totals = map { 0 } (0 .. scalar(@$alt));
 foreach my $sample (keys %$gt) {
   my $counts = $gt->{$sample}->{AD} or next;
   my @counts = split(/,/,$counts);
   for (my $i = 0; $i < scalar(@counts); $i++) { 
      $allele_totals[$i] += $counts[$i];
      $total += $counts[$i];
   } 
 }
 # Remove totally missing variants.
 @allele_totals = grep { $_ } @allele_totals;
 # Pass non-biallelic variants, we do not care about them.
 return 1 unless scalar(@allele_totals) != 2;
 my $thr = $self->thr * $total;
 return scalar(grep { $_ >= $thr } @allele_totals) == 2;
}



1;
