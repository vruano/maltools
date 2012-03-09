#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $out_file) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file ,
				'out=s' => \$out_file
) ;

my $usage = "USAGE : ./stats_typable.pl --manifest MANIFEST_FILE --vcf FILE --out FILE\n" ;
die $usage unless defined $manifest_file ;
die $usage unless defined $vcf_file ;
die $usage unless defined $out_file ;

read_manifest ( $manifest_file ) ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();
print $vcf->format_header() ;

open OUT , ">$out_file" ;
print OUT "#Chr\tPos\tRef\tAlt\tTypable\n" ;

while (my $x=$vcf->next_data_array()) {
	my $credible = 0 ;
	my $typable = 0 ;
	
	$typable = 1 if $x->[6] eq 'PASS' ; # Typable
	$typable = 1 if $x->[7] =~ m/\bLIT=/ ; # Literature SNP, counts as typable
	$credible = 1 if $x->[6] eq 'MinSamplesPerSNP' ; # Credible
	
	if ( $credible or $typable ) {
		my $o = $x->[0] . "\t" . $x->[1] . "\t" . $x->[3] . "\t" . $x->[4] . "\t" ;
		$o .= $typable ? 'Y' : 'N' ;
		print OUT "$o\n" ;
	}
	
	print join ( "\t" , @{$x} ) . "\n" ;
}

close OUT ;

$vcf->close() ;

0 ;
