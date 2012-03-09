#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my ( $manifest_file , $vcf_file , $range ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file ,
				'range=i' => \$range
) ;

my $usage = "USAGE : ./filter_coverage.pl --manifest FILE --vcf FILE --range NUMBER\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;
die $usage unless defined $range  ;

read_manifest ( $manifest_file ) ;

# Read coverage JSON object
my $cov ;
die "Coverage cutoff file does not exist\n" unless -e get_manifest_path ('coverage_cutoffs') ;
open COV , get_manifest_path ('coverage_cutoffs') ;
while ( <COV> ) {
	$cov .= $_ ;
}
close COV ;
$cov = decode_json ( $cov ) ;
die "Coverage cutoff file does not include range $range\n" unless defined $cov->{'ranges'}->{$range} ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'FILTER', ID=>'MinCoverage',Description=>"Low coverage, lower ${range}th centile (min ".$cov->{'ranges'}->{$range}->{'coding'}->[0]."x for coding, ".$cov->{'ranges'}->{$range}->{'noncoding'}->[0]."x for noncoding)" } ) ;
$vcf->add_header_line ( {key=>'FILTER', ID=>'MaxCoverage',Description=>"High coverage, upper ${range}th centile (max ".$cov->{'ranges'}->{$range}->{'coding'}->[2]."x for coding, ".$cov->{'ranges'}->{$range}->{'noncoding'}->[2]."x for noncoding)" } ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {
	my $c = $cov->{'positions'}->{$x->[0]}->[$x->[1]] || 0 ;
	

	if ( $x->[7] =~ m/\bCODING\b/ ) {
		set_filter_array ( $x , 'MinCoverage' ) if $c < $cov->{'ranges'}->{$range}->{'coding'}->[0] ;
		set_filter_array ( $x , 'MaxCoverage' ) if $c > $cov->{'ranges'}->{$range}->{'coding'}->[2] ;
	} else {
		set_filter_array ( $x , 'MinCoverage' ) if $c < $cov->{'ranges'}->{$range}->{'noncoding'}->[0] ;
		set_filter_array ( $x , 'MaxCoverage' ) if $c > $cov->{'ranges'}->{$range}->{'noncoding'}->[2] ;
	}
	
	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
