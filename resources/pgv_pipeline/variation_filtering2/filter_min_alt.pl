#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $p_total , $d_single ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'p_total=i' => \$p_total , 
				'd_single=i' => \$d_single
) ;

my $usage = "USAGE : ./filter_het_uniq.pl --manifest FILE --vcf FILE --p_total=MIN_PERCENT --d_single=MIN_DEPTH\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;
die $usage unless defined $p_total ;
die $usage unless defined $d_single ;

#read_manifest ( $manifest_file ) ;


my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'FILTER', ID=>'MinAlt',Description=>"SNP fails if the minor allele has not at least $p_total% in all samples, or more than ${d_single}x depth in at least one sample" } ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {

	die "No DP\n" unless $x->[7] =~ m/\bDP=(\d+)\b/ ;
	my $all_dp = $1 ;
	die "No AD\n" unless $x->[7] =~ m/\bAD=([0-9,]+)\b/ ;
	my @all_ad = split ',' , $1 ;
	
	my $minp_count = 0 ;
	if ( $all_dp > 0 ) {
		foreach ( @all_ad ) {
			$minp_count++ if $_ * 100 / $all_dp >= $p_total ;
		}
	}
	
	if ( $minp_count <= 1 ) { # None or only one allele >= x%, try per-sample absolute depth
		# Parse FORMAT field
		my @format = split ':' , $x->[8] ;
		my %format ;
		$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
		my $ad = $format{'AD'} ;
		die "No AD format set\n" unless defined $ad ;
	
		my $found_one = 0 ;
		my $max = scalar(@{$x})-1 ;
		my @min10 = ( 0 , 0 , 0 , 0 ) ;
		foreach ( 9 .. $max ) {
			my @d = split ':' , $x->[$_] ;
			my @calls = split ',' , $d[$ad] ;
#			my $cnt = 0 ;
			foreach my $n ( 0 .. $#calls ) {
				next unless $calls[$n] >= $d_single ;
				$min10[$n]++ ;
				$found_one = 1 if ($min10[0]>0?1:0) + ($min10[1]>0?1:0) + ($min10[2]>0?1:0) + ($min10[3]>0?1:0) > 1 ; # CAVEAT : Limited to 4 alleles. Might be a problem come indels.
#				$cnt++ if $c >= $d_single ;
			}
			last if $found_one ;
#			if ( $cnt > 1 ) {
#				$found_one = 1 ;
#				last ;
#			}
		}
		
		set_filter_array ( $x , 'MinAlt' ) if $found_one == 0 ;
	}

	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
