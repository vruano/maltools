#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my ( $manifest_file , $vcf_file , $use_coding , $use_noncoding ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'use_coding=i' => \$use_coding ,
				'use_noncoding=i' => \$use_noncoding
) ;

my $usage = "USAGE : ./filter_coding_and_regions.pl --manifest FILE --vcf FILE --use_coding=[0|1] --use_noncoding=[0|1]\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;

read_manifest ( $manifest_file ) ;

my %seq_ok ;
foreach ( @{$manifest->{'sequences'} || []} ) {
	$seq_ok{$_} = 1 ;
}

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

my $fd = [] ;
push @{$fd} , "coding" if $use_coding ;
push @{$fd} , "noncoding" if $use_noncoding ;
$fd = join ' and ' , @{$fd} ;
$vcf->add_header_line ( {key=>'FILTER', ID=>'CodingType',Description=>"Using $fd SNPs" } ) ;

$fd = join ',' , sort @{$manifest->{'sequences'}} ;
$vcf->add_header_line ( {key=>'FILTER', ID=>'Region',Description=>"Using SNPs from $fd" } ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {
	set_filter_array ( $x , 'Region' ) unless defined $seq_ok{$x->[0]} ;

	if ( $x->[7] =~ m/\bCODING\b/ ) {
		set_filter_array ( $x , 'CodingType' ) unless $use_coding ;
	} else {
		set_filter_array ( $x , 'CodingType' ) unless $use_noncoding ;
	}
	
	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
