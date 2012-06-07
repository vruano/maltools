#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my ( $manifest_file , $vcf_file , $min_cov ) ;
my $filter_name = 'MonoAllelic' ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'min_cov=i' => \$min_cov ,
				'filtername=s' => \$filter_name
) ;

my $usage = "USAGE : ./filter_monoallelic.pl --manifest FILE --vcf FILE --min_cov=MIN_COVERAGE [--filtername NAME|MonoAllelic]\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;
die $usage unless defined $min_cov ;

read_manifest ( $manifest_file ) ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'FILTER', ID=>$filter_name,Description=>"Fails if little evidence for alternate allele; min coverage=$min_cov" } ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {
	
	if ( $x->[7] =~ m/\bAD=\d+,\d+,\d+/ ) { # Tri-/quad-allelic state counts as evidence...
		# PASS
	} elsif ( $x->[7] =~ m/\bDP=0\b/ ) { # Blank SNP, no coverage
		set_filter_array ( $x , $filter_name ) ;
	} elsif ( $x->[7] !~ m/\bAD=\d+,/ ) { # Mono-allelic reference
		set_filter_array ( $x , $filter_name ) ;
	} elsif ( $x->[7] =~ m/\bAD=(\d+),(\d+)/ ) { # Bi-allelic
		
		my ( $ref , $alt ) = ( $1 , $2 ) ;
		if ( $ref < $min_cov ) { # Almost all alt
			set_filter_array ( $x , $filter_name ) ;
		} elsif ( $alt < $min_cov ) { # Almost all ref
			set_filter_array ( $x , $filter_name ) ;
		} else {
			
			# Now for some real work...
			my @format = split ':' , $x->[8] ;
			my %format ;
			$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
			my $ad = $format{'AD'} ;
			die "No AD format set\n" unless defined $ad ;
			
			my @cnt = ( 0 , 0 ) ;
			my $max = scalar(@{$x})-1 ;
			foreach ( 9 .. $max ) {
				my @d = split ':' , $x->[$_] ;
				die "Bad format\n" unless $d[$ad] =~ m/^(\d+),(\d+)$/ ;
				$cnt[0]++ if $1 >= $min_cov ;
				$cnt[1]++ if $2 >= $min_cov ;
				last if $cnt[0] > 0 and $cnt[1] > 0 ;
			}
			
			unless ( $cnt[0] > 0 and $cnt[1] > 0 ) {
				set_filter_array ( $x , $filter_name ) ;
			}
			
		}
	} # PASS by default

	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
