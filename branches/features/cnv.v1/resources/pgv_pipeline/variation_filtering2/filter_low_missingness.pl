#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $min_snps , $min_samples , $json_file , $min_cov ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'min_snps=i' => \$min_snps , 
				'min_samples=i' => \$min_samples , 
				'json=s' => \$json_file ,
				'min_cov=i' => \$min_cov
) ;

my $usage = "USAGE : ./filter_low_missingness.pl --manifest FILE --vcf FILE --json FILE --min_snps MIN_SNPS_PER_SAMPLE --min_samples MIN_SAMPLES_PER_SNP --min_cov MIN_COVERAGE\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;
die $usage unless defined $json_file ;
die $usage unless defined $min_snps ;
die $usage unless defined $min_samples ;
die $usage unless defined $min_cov ;

read_manifest ( $manifest_file ) ;

my $data = read_json_from_file ( $json_file ) ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

my (@samples) = $vcf->get_samples();

my @keep_samples ;
my @keep_samples_id ;
foreach my $sid ( 0 .. $#samples ) {
	my $sample = $samples[$sid] ;
	next unless defined $data->{'snps_in_samples'}->{$sample} ;
	my $sum = 0 ;
	$sum += $_ foreach values %{$data->{'snps_in_samples'}->{$sample}} ;
	next if $sum < $min_snps ;
	push @keep_samples , $sample ;
	push @keep_samples_id , $sid ;
}

my $samples_removed = scalar ( @samples ) - scalar ( @keep_samples ) ;
$vcf->add_header_line ( {key=>'FILTER', ID=>'MinSamplesPerSNP',Description=>"SNP fails if it shows < $min_samples times at ${min_cov}x" } ) ;
$vcf->add_header_line ( {key=>'FILTER', ID=>'MinSNPsPerSample',Description=>"Removed $samples_removed samples with <$min_snps at ${min_cov}x" } ) ;
print $vcf->format_header(\@keep_samples);

# Parse through data
while (my $x=$vcf->next_data_array()) {

	my @o ;
	push @o , $x->[$_] foreach ( 0 .. 8 ) ;
	push @o , $x->[9+$_] foreach @keep_samples_id ;

	my $id = $o[0] . "\t" . $o[1] ;
	
	my @format = split ':' , $o[8] ;
	my %format ;
	$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
	my $dp = $format{'DP'} ;
	die "No DP format set\n" unless defined $dp ; # Paranoia
	
	my $sample_count = 0 ;
#	foreach my $sid ( 9 .. $#o ) {
#		my @parts = split ':' , $o[$sid] ;
	my $last_x = scalar ( @{$x} ) - 1 ;
	foreach my $sid ( 9 .. $last_x ) {
		my @parts = split ':' , $x->[$sid] ;
		my $depth = $parts[$dp] ;
		$sample_count++ if $depth >= $min_cov ;
	}
	
	set_filter_array ( \@o , 'MinSamplesPerSNP' ) if $sample_count < $min_samples ;

	$o[6] = 'PASS' if $o[6] eq '.' ;
	
	
	print join ( "\t" , @o ) . "\n" ;
}

$vcf->close() ;

0 ;
