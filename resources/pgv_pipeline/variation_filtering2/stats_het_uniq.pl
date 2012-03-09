#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my ( $manifest_file , $vcf_file , $min_cov ) ;

GetOptions  (
				'vcf=s' => \$vcf_file , 
) ;

my $usage = "USAGE : ./filter_biallelic.pl --manifest FILE --vcf FILE --min_cov=MIN_COVERAGE\n" ;
die $usage unless defined $vcf_file ;

#read_manifest ( $manifest_file ) ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

# Parse through data
my @uniq ;
while (my $x=$vcf->next_data_array()) {

#	next unless $x->[6] eq 'PASS' ; # ???
	
	next unless $x->[7] =~ m/\bAC=([0-9,]+)\b/ ;
	my @d = split ',' , $1 ;
	next unless 1 < scalar @d ;
	
	next unless $x->[7] =~ m/\bUQ=(\d+)\b/ ;
	my $uq = $1 ;

	# Parse FORMAT field
	my @format = split ':' , $x->[8] ;
	my %format ;
	$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
	my $ad = $format{'AD'} ;
	die "No AD format set\n" unless defined $ad ;

	my $het_samples = 0 ;
	my $max = scalar(@{$x})-1 ;
	my $nos = 0 ;
	foreach ( 9 .. $max ) {
		my @d = split ':' , $x->[$_] ;
		my @calls = split ',' , $d[$ad] ;
		my $cnt = 0 ;
		foreach my $c ( @calls ) {
			$cnt++ if $c > 0 ;
		}
		if ( $cnt > 0 ) {
			$nos++ ;
			$het_samples++ if $cnt > 1 ;
		}
	}
	
	push @{$uniq[$uq]} , $het_samples / $nos if $nos > 0 ;

}

$vcf->close() ;

print "Uniqueness\t#het\n" ;
foreach ( 0 .. $#uniq ) {
	my @a = sort { $a <=> $b } @{$uniq[$_] || [] } ;
	if ( 0 == scalar @a ) {
		print "$_\t0\n" ;
	} else {
		print "$_\t" . ($a[int($#a/2)]) . "\n" ;
	}
}

0 ;
