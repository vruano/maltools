#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $snp_list ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'snps=s' => \$snp_list
) ;

my $usage = "USAGE : ./filter_pass_literature_snps.pl --manifest FILE --vcf FILE [--snps SNP_LIST_FILE]\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;

read_manifest ( $manifest_file ) ;

# SNP list
$snp_list = get_manifest_path ( "literature_snps" ) unless defined $snp_list ;
die $usage unless defined $snp_list ;
die "SNP file does not exist : $snp_list\n" unless -e $snp_list ;

my %ls ;
open FILE , $snp_list ;
while ( <FILE> ) {
	chomp ;
	next unless $_ =~ m/^(\S+)\s(\d+)/ ;
	$ls{$1}->{$2} = 1 ;
}
close FILE ;

# VCF header
my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'INFO', ID=>'LIT',Number=>'1',Type=>'String',Description=>'Literature SNP'} ) ;

print $vcf->format_header();


# Parse through data
while (my $x=$vcf->next_data_array()) {

	if ( defined $ls{$x->[0]}->{$x->[1]} ) {
		my $comment = "." ; # Comment about the literature SNP, for later...
		$x->[7] .= ';LIT=' . $comment ;
	}

	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
