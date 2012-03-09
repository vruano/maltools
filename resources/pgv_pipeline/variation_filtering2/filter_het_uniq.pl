#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $min_cov , $max_uniq ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'min_cov=i' => \$min_cov , 
				'max_uniq=i' => \$max_uniq
) ;

my $usage = "USAGE : ./filter_het_uniq.pl --manifest FILE --vcf FILE --min_cov=MIN_COVERAGE --max_uniq=MAX_UNIQUENESS_FOR_HET\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;
die $usage unless defined $min_cov ;
die $usage unless defined $max_uniq ;

read_manifest ( $manifest_file ) ;


my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'FILTER', ID=>'HetUniq',Description=>"SNP fails if it has heterozygous samples(min coverage=$min_cov) and a uniqueness score>$max_uniq" } ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {

	$x->[7] =~ m/\bUQ=(\d+)\b/ ;
	my $uniqueness = $1 ;
	
	if ( $x->[7] =~ m/\bAD=\d+,\d+/ and $uniqueness > $max_uniq ) { # Not mono-allelic
	
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
			if ( $cnt > 1 ) {
				$nos++ ;
				$het_samples++ if $cnt >= $min_cov ;
				last if $het_samples > 0 ;
			}
		}
		
		set_filter_array ( $x , 'HetUniq' ) if $het_samples > 0 ;

	}

	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
