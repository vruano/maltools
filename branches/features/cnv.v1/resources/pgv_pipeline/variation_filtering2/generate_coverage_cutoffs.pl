#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;

my $production = 0 ;

my ( $manifest_file , $data_file ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'datafile=s' => \$data_file # To reuse an existing JSON file, but recalculate the regions
) ;

die "Needs --manifest\n" unless defined $manifest_file  ;
read_manifest ( $manifest_file ) ;

my %seq_ok ;
foreach ( @{$manifest->{'sequences'} || []} ) {
	$seq_ok{$_} = 1 ;
}

# Coding regions
my %is_coding ;
print STDERR "Reading coding regions...\n" unless $production ;
open FILE , get_manifest_path('coding_regions') ;
while ( <FILE> ) {
	next unless $_ =~ m/^(\S+)\s(\d+)\s(\d+)$/ ;
	my ( $chr , $from , $to ) = ( $1 , $2 , $3 ) ;
	next unless defined $seq_ok{$chr} ;
	foreach my $p ( $from .. $to ) {
		$is_coding{$chr}->[$p] = 1 ;
	}
}
close FILE ;

my %results ;
$results{'samples'} = [] ;
$results{'positions'} = {} ;

# Total coverage
print STDERR "Reading coverage data...\n" unless $production ;
if ( defined $data_file ) {
	my $cov ;
	open COV , get_manifest_path ('coverage_cutoffs') ;
	while ( <COV> ) {
		$cov .= $_ ;
	}
	close COV ;
	%results = %{decode_json ( $cov )} ;
	delete $results{'ranges'} if defined $results{'ranges'} ;
} else {
	foreach my $sample ( @{$manifest->{'samples'}} ) {
		my $file = $manifest->{'paths'}->{'coverage_path'} . '/' . $sample . '.coverage' ;
		unless ( -e $file ) {
                       die "ERROR $sample - coverage file not found\n" ;
			next ;
		}
		print "Processing $sample\n" ;
		push @{$results{'samples'}} , $sample ;
		open FILE , "gawk '{ print \$1 \"\\t\" \$2 \"\\t\" (\$4+\$5+\$6+\$7) }' $file |" ;
		my $rp = $results{'positions'} ;
		while ( <FILE> ) {
			next unless $_ =~ m/^(\S+)\s(\d+)\s(\d+)$/ ;
			next unless defined $is_coding{$1} ;
			$rp->{$1}->[$2] += $3 ;
		}
		close FILE ;
	}
}


# Count coverage occurrences
print STDERR "Counting occurrences...\n" unless $production ;
my @coding ;
my @noncoding ;
foreach my $chr ( keys %{$results{'positions'}} ) {
	my $max = scalar @{$results{'positions'}->{$chr}} ;
	foreach my $pos ( 1 .. $max ) {
		if ( defined $is_coding{$chr}->[$pos] ) {
			$coding[$results{'positions'}->{$chr}->[$pos]||0]++ ;
		} else {
			$noncoding[$results{'positions'}->{$chr}->[$pos]||0]++ ;
		}
	}
}


# Calculate ranges
foreach ( 0 .. 9 ) {
	my %r ;
	my $range = 5 * $_ ;
	$r{'range'} = [ $range , 100 - $range ] ;
	$r{'coding'} = get_range ( \@coding , $range ) ;
	$r{'noncoding'} = get_range ( \@noncoding , $range ) ;
	print STDERR encode_json ( \%r ) . "\n" ; # DEBUG
	$results{'ranges'}->{$range} = \%r ;
}

open OUT , ">" . get_manifest_path('coverage_cutoffs') ;
print OUT encode_json ( \%results ) ;
close OUT ;

0 ;

sub get_range {
	my ( $a , $range ) = @_ ;
	my @ret = ( 0 , 0 , 0 ) ;
	my $sum = 0 ;
	$sum += ( $_ || 0 ) foreach @{$a} ;
	
	my @bound = ( $sum*$range/100 , $sum/2 , $sum-$sum*$range/100 ) ;
	my $cnt = 0 ;
	my $max = scalar ( @{$a} ) - 1 ;
	foreach my $p ( 0 .. $max ) {
		my $cnt_new = $cnt + ( $a->[$p] || 0 ) ;
		$ret[0] = $p if $cnt <= $bound[0] and $cnt_new > $bound[0] ;
		$ret[1] = $p if $cnt <= $bound[1] and $cnt_new > $bound[1] ;
		$ret[2] = $p if $cnt <= $bound[2] and $cnt_new > $bound[2] ;
		$cnt = $cnt_new ;
	}
	return \@ret ;
}
