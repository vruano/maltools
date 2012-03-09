package MalariaGEN::AGV::Vcf::Filters::BIAL;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::Vcf::Filter';

has '+id' => ( default => 'BIAL' );
has '+description' => ( default => 'Biallelic filter');

sub passes {
 my ($self,$h) = @_; 

 my $alt = $h->{ALT};
 my $gt = $h->{gtypes};
 return 0 if scalar(@$alt) > 2;
 my @allele_totals = map { 0 } (0 .. scalar(@$alt));
 foreach my $sample (keys %$gt) {
   my $counts = $gt->{$sample}->{AD} or next;
   my @counts = split(/,/,$counts);
   for (my $i = 0; $i < scalar(@counts); $i++) { 
      $allele_totals[$i] += $counts[$i];
   } 
 }
 return 0 if 2 != scalar(grep { $_ > 1 } @allele_totals);
 return 1;
}



1;
