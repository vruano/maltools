#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my ( $manifest_file , $vcf_file , $file_base , $steps ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'steps=s' => \$steps , 
				'base=s' => \$file_base
) ;

my $usage = "USAGE : ./stats_low_missingness.pl --vcf FILE --base OUT_FILE_BASENAME --steps 5,7,10,15,...\n" ;
die $usage unless defined $vcf_file ;
die $usage unless defined $file_base ;
die $usage unless defined $steps ;

my @steps = split ',' , $steps ;

#read_manifest ( $manifest_file ) ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();
print $vcf->format_header();

my (@samples) = $vcf->get_samples();
my %snps_in_samples ;
my %samples_per_snp ;

# Parse through data
while (my $x=$vcf->next_data_array()) {
	print join ( "\t" , @{$x} ) . "\n" ;

	next unless $x->[6] eq 'PASS' ; # Survivors only
	
	my $id = $x->[0] . "\t" . $x->[1] ;
	
	my @format = split ':' , $x->[8] ;
	my %format ;
	$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
	my $dp = $format{'DP'} ;
	die "No DP format set\n" unless defined $dp ;
	
	foreach my $sid ( 0 .. $#samples ) {
		my @parts = split ':' , $x->[9+$sid] ;
		my $depth = $parts[$dp] ;
		my $max = 0 ;
		foreach ( @steps ) {
			if ( $depth >= $_ ) {
				$max = $_ ;
			} else {
				last ;
			}
		}
		
		if ( $max > 0 ) {
			$snps_in_samples{$samples[$sid]}->{$max}++ ;
			$samples_per_snp{$id}->{$max}++ ;
		}
		
	}

}

$vcf->close() ;

# JSON
my %out ;
$out{'snps_in_samples'} = \%snps_in_samples ;
#$out{'samples_per_snp'} = \%samples_per_snp ;
open FILE , ">$file_base.json" ;
print FILE encode_json ( \%out ) ;
close FILE ;

write_sorted ( "$file_base.snps_in_samples.tab" , \%snps_in_samples , 1000 ) ;
write_sorted ( "$file_base.samples_per_snp.tab" , \%samples_per_snp , 1000 ) ;

0 ;

sub write_sorted {
	my ( $filename , $data , $max ) = @_ ;
	my @d ;
	foreach my $v ( values %{$data} ) {
		my @n = ( 0 ) ;
		foreach ( @steps ) {
			push @n , ( $v->{$_} || 0 ) ;
			$n[0] += ( $v->{$_} || 0 ) ;
		}
		push @d , \@n ;
	}
	@d = sort { $a->[0] <=> $b->[0] } @d ;
	
	# Subset if too large
	if ( $max < scalar @d ) {
		my $diff = int ( scalar ( @d ) / $max ) + 1 ;
		my @d2 ;
		foreach ( 0 .. $#d ) {
			push @d2 , $d[$_] if $_ % $diff == 0 ;
		}
		@d = @d2 ;
	}
	
	# Write to tabbed file
	open FILE , ">$filename" ;
	print FILE '#' . join ( "\t#" , @steps ) . "\tSum\n" ;
	foreach my $x ( @d ) {
		push @{$x} , shift @{$x} ;
		print FILE join ( "\t" , @{$x} ) . "\n" ;
	}
	close FILE ;
}
