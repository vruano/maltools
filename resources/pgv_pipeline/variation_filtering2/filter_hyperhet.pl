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

my $usage = "USAGE : ./filter_hyperhet.pl --manifest FILE --vcf FILE\n" ;
die $usage unless defined $manifest_file  ;
die $usage unless defined $vcf_file ;

read_manifest ( $manifest_file ) ;

# Group data
my @groups = sort keys %{$manifest->{'groups'}} ;
my %is_hh ;
my $hh_dir = $manifest->{'paths'}->{'out_path'} . '/hyperheterozygous' ;
foreach my $group ( @groups ) {
	my $file = "$hh_dir/hyperhet.$group.pls" ;
	die "PLS file for group $group not found : $file\n" unless -e $file ;
	open FILE , "cut -f 1,2,16 $file |" ;
	while ( <FILE> ) {
		chomp ;
		my ( $chr , $pos , $pls ) = split "\t" , $_ ;
		next unless defined $pls ;
		next unless $pls =~ m/^[0-9.]+$/ ;
		next if $pls <= $manifest->{'groups'}->{$group}->{'hyperhet_cutoff'} ;
		push @{$is_hh{$chr}->{$pos}} , $group ; # Maybe add groups to info tag later...
	}
	close FILE ;
}

# VCF header
my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

$vcf->add_header_line ( {key=>'FILTER', ID=>'HyperHet',Description=>"SNP fails if it is hyper-heterozygous in one of the groups (" . join(',',@groups) . "); filter is applied to low-missingness SNPs only" } ) ;

print $vcf->format_header();

# Parse through data
while (my $x=$vcf->next_data_array()) {
	set_filter_array ( $x , 'HyperHet' ) if defined $is_hh{$x->[0]}->{$x->[1]} ;
	$x->[6] = 'PASS' if $x->[6] eq '.' ;
	print join ( "\t" , @{$x} ) . "\n" ;
}

$vcf->close() ;

0 ;
