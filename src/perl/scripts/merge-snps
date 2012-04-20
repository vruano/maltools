#!/usr/bin/perl

use strict ;
use warnings ;

my ( %ref , %alt ) ;
my $read = 0 ;
my %iupac = (
	'A' => 'A',
	'C' => 'C',
	'G' => 'G',
	'T' => 'T',
	'R' => 'AG',
	'Y' => 'CT',
	'M' => 'AC',
	'K' => 'GT',
	'W' => 'AT',
	'S' => 'CG',
	'B' => 'CGT',
	'D' => 'AGT',
	'H' => 'ACT',
	'V' => 'ACG',
	'N' => 'ACGT'
) ;
my %iupac_reverse ;
foreach ( keys %iupac ) {
	my $v = $iupac{$_} ;
	$iupac_reverse{$v} = $_ ;
}

while ( <> ) {
	$_ =~ /^([\S]+)\s+(\d+)\s+(.)\s+(.)\s*$/ ;
	next unless defined $4 ;
	$read++ ;
	
	print STDERR "REFERENCE MISMATCH!\n" if ( $ref{$1}->[$2] || $3 ) ne $3 ;
	
	$ref{$1}->[$2] = $3 ;
	$alt{$1}->[$2] = merge_codes ( $4 , $alt{$1}->[$2] || $4 ) ;
}

my $out = 0 ;
foreach my $chr ( sort keys %alt ) {
	my @data = @{$alt{$chr}} ;
	foreach ( 0 .. $#data ) {
		next unless defined $data[$_] ;
		print "$chr\t$_\t" . $ref{$chr}->[$_] . "\t" . $alt{$chr}->[$_] . "\n" ;
		$out++ ;
	}
}

print STDERR "Read $read SNPs\n" ;
print STDERR "Wrote $out SNPs\n" ;

sub merge_codes {
	my ( $n , $o ) = @_ ;
	return $n if $n eq $o ;
	my $new_all = $iupac{$n} . $iupac{$o} ;
	my $k = '' ;
	$k .= 'A' if $new_all =~ m/A/ ;
	$k .= 'C' if $new_all =~ m/C/ ;
	$k .= 'G' if $new_all =~ m/G/ ;
	$k .= 'T' if $new_all =~ m/T/ ;
	
#	print "$n\t$o\t$new_all\t$k\t" . $iupac_reverse{$k} . "\n" if $iupac_reverse{$k} ne $n ;
	return $iupac_reverse{$k} ;
}
