#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $min_cov ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'min_cov=i' => \$min_cov
) ;

my $usage = "USAGE : ./filter_biallelic.pl --manifest FILE --vcf FILE --min_cov=MIN_COVERAGE\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;
die $usage unless defined $min_cov ;

read_manifest ( $manifest_file ) ;


my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'FILTER', ID=>'Biallelic',Description=>"Only biallelic SNPs allowed; low-coverage tri-allelic SNPS may be converted to biallelic, if the allele coverage is less than $min_cov" } ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {
	
	if ( $x->[7] =~ m/\bAD=\d+,\d+,\d+/ ) { # Tri-/quad-allelic state
	
		$x->[7] =~ m/\bAD=([0-9,]+)\b/ ;
		my @nums = split ',' , $1 ;
		
		if ( 3 != scalar @nums ) { # Not triallelic
			set_filter_array ( $x , 'Biallelic' ) ;
		} else {
			my @remove ;
			foreach ( 0 .. $#nums ) {
				push @remove , $_ if $nums[$_] < $min_cov ;
			}
			if ( 1 != scalar @remove ) { # Not a single clear candidate to remove
				set_filter_array ( $x , 'Biallelic' ) ;
			} elsif ( $remove[0] == 0 ) { # Remove reference == bad
				set_filter_array ( $x , 'Biallelic' ) ;
			} else {
				splice @nums , $remove[0] , 1 ;
				
				# Fix ALT field
				my @a = split ',' , $x->[4] ;
				splice @a , $remove[0]-1 , 1 ;
				$x->[4] = join ',' , @a ;
				
				# Fix INFO field
#				print $x->[7] . " => " ;
				my $j = join ',' , @nums ;
				$x->[7] =~ s/\bAD=[0-9,]+\b/AD=$j/ ;
				$j = 0 ;
				$j += $_ foreach @nums ;
				$x->[7] =~ s/\bDP=\d+\b/DP=$j/ ;
#				print $x->[7] . "\n" ;
				
				# Parse FORMAT field
				my @format = split ':' , $x->[8] ;
				my %format ;
				$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
				my $ad = $format{'AD'} ;
				my $dp = $format{'DP'} ;
				
				my $max = scalar(@{$x})-1 ;
				foreach my $s ( 9 .. $max ) {
#					print $x->[$s]  . ' => ' ;
					my @d = split ':' , $x->[$s] ;
					my @calls = split ',' , $d[$ad] ;
					splice @calls , $remove[0] , 1 ;
					my $sum = 0 ;
					$sum += $_ foreach @calls ;
					$x->[$s] =~ s/\bDP=\d+\b/DP=$sum/ ;
					$d[$ad] = join ',' , @calls ;
					$x->[$s] = join ':' , @d ;
#					print $x->[$s]  . '; ' ;
				}
#				print "\n" ;

				# TODO fix DP
				
				# Annotate INFO field
				@a = split ';' , $x->[7] ;
				push @a , 'TRI2BI' ;
				$x->[7] = join ';' , @a ;
			}
		}
		
		print join ( "\t" , @{$x} ) . "\n" if  $x->[6] =~ m/Biallelic/ and not $production ;
	}

	next unless $production ;

	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
