#!/software/bin/perl

#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)
# Group:      Team 112
# SVN:        malaria, path /pipelines/magnus
#

use strict ;
use warnings ;
use Getopt::Long ;
use JSON::XS ;
use Data::Dumper ;

package Variation ;

use strict ;
use warnings ;
use Getopt::Long ;
use JSON ;
use Data::Dumper ;
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(
    $manifest
    $outdir

	%iupac
	%iupac_reverse
	%amino_acid
	%codon2aa

	read_json_from_file
    read_manifest
    get_manifest_path
    initialize_iupac
    set_filter_hash
    set_filter_array
) ;

our $manifest ;
our $outdir ;

our %iupac ;
our %iupac_reverse ;
our %amino_acid ;
our %codon2aa ;

# Initialize
initialize_iupac() ;

sub read_json_from_file {
    my ( $file ) = @_ ;
    my $s ;
    open FILE , $file ;
    $s .= $_ foreach ( <FILE> ) ;
    close FILE ;
    return from_json ( $s ) ;
}

sub read_manifest {
    my ( $file ) = @_ ;
    $manifest = read_json_from_file ( $file ) ;
    
    $outdir = $manifest->{'paths'}->{'out_path'} ;
    die "No output path defined\n" unless defined $outdir ;
    mkdir $outdir unless -d $outdir ;
}

sub get_manifest_path {
    my ( $k ) = @_ ;
    my $p = $manifest->{'paths'}->{$k} || die "Path to $k not defined in manifest\n" ;
    return $p if $p =~ m|^[/.]| ;
    return $manifest->{'paths'}->{'base_path'} . '/' . $p ;
}

sub initialize_iupac {
	# Initialize IUPAC codes
	$iupac{'A'} = 'A' ;
	$iupac{'C'} = 'C' ;
	$iupac{'G'} = 'G' ;
	$iupac{'T'} = 'T' ;
	$iupac{'AG'} = 'R' ;
	$iupac{'CT'} = 'Y' ;
	$iupac{'AC'} = 'M' ;
	$iupac{'GT'} = 'K' ;
	$iupac{'AT'} = 'W' ;
	$iupac{'CG'} = 'S' ;
	$iupac{'CGT'} = 'B' ;
	$iupac{'AGT'} = 'D' ;
	$iupac{'ACT'} = 'H' ;
	$iupac{'ACG'} = 'V' ;
	$iupac{'ACGT'} = 'N' ;

	foreach ( keys %iupac ) {
		$iupac_reverse{$iupac{$_}} = $_ ;
	}
	
	# Initialize amino acids
	$amino_acid{'Ala/A'} = [ 'GCT', 'GCC', 'GCA', 'GCG' ] ;
	$amino_acid{'Leu/L'} = [ 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG' ] ;
	$amino_acid{'Arg/R'} = [ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ] ;
	$amino_acid{'Lys/K'} = [ 'AAA', 'AAG' ] ;
	$amino_acid{'Asn/N'} = [ 'AAT', 'AAC' ] ;
	$amino_acid{'Met/M'} = [ 'ATG' ] ;
	$amino_acid{'Asp/D'} = [ 'GAT', 'GAC' ] ;
	$amino_acid{'Phe/F'} = [ 'TTT', 'TTC' ] ;
	$amino_acid{'Cys/C'} = [ 'TGT', 'TGC' ] ;
	$amino_acid{'Pro/P'} = [ 'CCT', 'CCC', 'CCA', 'CCG' ] ;
	$amino_acid{'Gln/Q'} = [ 'CAA', 'CAG' ] ;
	$amino_acid{'Ser/S'} = [ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ] ;
	$amino_acid{'Glu/E'} = [ 'GAA', 'GAG' ] ;
	$amino_acid{'Thr/T'} = [ 'ACT', 'ACC', 'ACA', 'ACG' ] ;
	$amino_acid{'Gly/G'} = [ 'GGT', 'GGC', 'GGA', 'GGG' ] ;
	$amino_acid{'Trp/W'} = [ 'TGG' ] ;
	$amino_acid{'His/H'} = [ 'CAT', 'CAC' ] ;
	$amino_acid{'Tyr/Y'} = [ 'TAT', 'TAC' ] ;
	$amino_acid{'Ile/I'} = [ 'ATT', 'ATC', 'ATA' ] ;
	$amino_acid{'Val/V'} = [ 'GTT', 'GTC', 'GTA', 'GTG' ] ;
	$amino_acid{'STP/X'} = [ 'TAG', 'TGA', 'TAA' ] ;

	foreach my $aa ( keys %amino_acid ) {
		my $aa2 = substr $aa , 4 , 1 ;
		foreach ( @{$amino_acid{$aa}} ) {
			$codon2aa{$_} = $aa2 ;
		}
	}
	
}

sub set_filter_hash {
	my ( $x , $filter ) = @_ ;
	$x->{'FILTER'}->[0] = 'PASS' if ( 1 == scalar ( @{$x->{'FILTER'}} ) and $x->{'FILTER'}->[0] eq '.' ) ;
	foreach ( @{$x->{'FILTER'}} ) {
		return if $_ eq $filter ;
	}
	if ( 1 == scalar ( @{$x->{'FILTER'}} ) and $x->{'FILTER'}->[0] eq 'PASS' ) {
		$x->{'FILTER'}->[0] = $filter ;
	} else {
		push @{$x->{'FILTER'}} , $filter ;
	}
}

sub set_filter_array {
	my ( $x , $filter ) = @_ ;
	return if $x->[6] =~ m/\b($filter)\b/ ;
	$x->[6] = 'PASS' if '.' eq $x->[6] ;
	if ( $x->[6] eq 'PASS' ) {
		$x->[6] = $filter ;
	} else {
		$x->[6] .= ';' . $filter ;
	}
}


1 ;
