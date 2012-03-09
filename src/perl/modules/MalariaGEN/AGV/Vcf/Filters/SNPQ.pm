package MalariaGEN::AGV::Vcf::Filters::SNPQ;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::Vcf::Filter';

has '+id' => ( default => 'SNPQ' );
has '+description' => ( default => 'SNP Quality filter Anopheles edition');

sub passes {
 my ($self,$h) = @_; 
 my $info = $h->{INFO};

 #  * the base pair at the position is not reference in all samples,
 # INTRINSIC no need to check.

 #   * the base pair at the position is not non-reference in all samples
 # DUE to bad reference quality and working across species barrier make no sense for Anopheles.
 # Refrased as this: at least there are two alternative alleles or the reference is present in
 # in more than one read and so on. Basically we allow for an alternative to take on the references role.
 my $alt = $h->{ALT};
 my $gt = $h->{gtypes};
 my @allele_totals = map { 0 } (0 .. scalar(@$alt));
 my @allele_maxs = @allele_totals;
 foreach my $sample (keys %$gt) {
   my $counts = $gt->{$sample}->{AD} or next;
   my @counts = split(/,/,$counts);
   for (my $i = 0; $i < scalar(@counts); $i++) { 
      $allele_totals[$i] += $counts[$i];
      $allele_maxs[$i] = $counts[$i] if $counts[$i] > $allele_maxs[$i];
   } 
 }
 #  * the position has the reference base for more than 1 read over all samples,
 #  * the position has the non-reference base for more than 1 read over all samples,
 return 0 unless  2 <= scalar(grep { $_ > 1 } @allele_totals);
 #   * the position has the reference base for more than 1 read within at least one sample,
 #   * and, the position has the non-reference base for more than 1 read within at least one sample, 
 return 0 unless 2 <= scalar(grep { $_ > 1} @allele_maxs);
 return 1;
}



1;
