#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my ( $manifest_file , $vcf_file , $samples , $keep ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'samples=s' => \$samples , 
				'keep=i' => \$keep
) ;

my $usage = "USAGE : ./remove_samples.pl --vcf FILE --samples \"SAMPLE1,SAMPLE2,...\" [--keep 1]\n" ;
die $usage unless defined $vcf_file ;
die $usage unless defined $samples ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

my (@current_samples) = $vcf->get_samples();

my %new_samples ;
$new_samples{$_} = 1 foreach ( split ',' , $samples ) ;

my @removed_samples ;
my @keep_samples ;
my @keep_samples_id ;
foreach my $sid ( 0 .. $#current_samples ) {
	my $sample = $current_samples[$sid] ;
	if ( ( $keep and defined $new_samples{$sample} ) or ( 1 != ($keep||0) and not defined $new_samples{$sample} ) ) {
		push @keep_samples , $sample ;
		push @keep_samples_id , $sid ;
	} else {
		push @removed_samples , $sample ;
	}
}

my $rs = join ',' , @removed_samples ;
$vcf->add_header_line ( {key=>'FILTER', ID=>'SampleRemove',Description=>"Removed samples $rs" } ) ;
print $vcf->format_header(\@keep_samples);

# Parse through data
while (my $x=$vcf->next_data_array()) {

	my @o ;
	push @o , $x->[$_] foreach ( 0 .. 8 ) ;
	push @o , $x->[9+$_] foreach @keep_samples_id ;
	
	print join ( "\t" , @o ) . "\n" ;
}

$vcf->close() ;

0 ;
