#!/usr/bin/perl

#########
# Author:     Magnus Manske (mm6@sanger.ac.uk)
# Group:      Team 112
# SVN:        malaria, path /pipelines/magnus
#

package gff2json ;

require Exporter;
use strict ;
use warnings ;
use URI::Escape;

our @ISA = qw(Exporter);
our @EXPORT = qw(
    get_gff_object_from_file
    get_gff_object_from_filehandle
    read_fasta
) ;

# To "unify" chromosome neames
my %hard_replace = ( 'NC_002375' => 'MITO1' , 'X95275' => 'PLAST1' , 'X95276' => 'PLAST2' , 
                    'API_IRAB' => '' , 'M76611' => '' , 'PfNF54' => '' , 'PFC10_API_IRAB' => '' ) ;

sub get_gff_object_from_file {
    my ( $file ) = @_ ;
    my $fh ;
    open $fh , $file ;
    my $ret = get_gff_object_from_filehandle ( $fh ) ;
    close $fh ;
    return $ret ;
}

sub get_gff_object_from_filehandle {
    my ( $fh ) = @_ ;
    my %json ;
    my %hadthat ;
    while ( <$fh> ) {
        next if $_ =~ m/^#/ ;
        chomp ;
        my ( $seqid , $source , $type , $from , $to , $score , $strand , $phase , $attr ) = split "\t" , $_ ;
        next unless defined $attr ;
        
        $seqid =~ s/^.+\|// ; # Remove "apidb|" prefix
        $type = lc $type ; # Paranoia
        
        my %d ;
        $d{'type'} = $type ;
        $d{'from'} = $from ;
        $d{'to'} = $to ;
        $d{'strand'} = $strand ;
        
        foreach my $a ( split ';' , $attr ) {
            next unless $a =~ m/^([^=]+)=(.+)$/ ;
            $d{lc $1} = $2 ;
        }
        foreach my $ot ( 'ontology_term' , 'dbxref' ) {
            next unless defined $d{$ot} ;
            my @x = split ',' , $d{$ot} ;
            $d{$ot} = \@x ;
        }
    
        next unless defined $d{'id'} ;
    
        foreach ( keys %d ) {
            next if $_ eq 'parent' ;
            next if ref($d{$_}) eq 'ARRAY' ;
            $d{$_} =~ tr/+/ / ;
            $d{$_} = uri_unescape ( $d{$_} ) ;
        }
        $d{'id'} =~ s/^.+\|// ; # Remove "apidb|" prefix
        $d{'parent'} =~ s/^.+\|// if defined $d{'parent'} ; # Remove "apidb|" prefix
        
        my $id = $d{'id'} ;
        delete $d{'id'} ;
        push @{$json{$seqid}->{$type}} , $id unless defined $hadthat{$seqid}->{$type}->{$id} ;
        $hadthat{$seqid}->{$type}->{$id} = 1 ;

        push @{$json{$seqid}->{'data'}->{$id}} , \%d ;
    }
    
    # Hackish fixes to unify annotation and reference
    $hard_replace{sprintf 'Pf3D7_%02d',$_} = 'MAL'.$_ foreach ( 1 .. 14 ) ;
    foreach ( keys %hard_replace ) {
        next unless defined $json{$_} ;
        $json{$hard_replace{$_}} = $json{$_} unless '' eq $hard_replace{$_} ;
        delete $json{$_} ;
    }
    
    return \%json ;
}

sub read_fasta {
	my ( $file ) = @_ ;
	my %ret ;
	my $seq ;
	my $name ;
	open FILE , $file ;
	while ( <FILE> ) {
		chomp ;
		if ( $_ =~ /^>/ ) {
			if ( defined $seq ) {
			    $seq =~ s/\s//g ;
				$ret{$name} = $seq ;
			}
			$name = substr $_ , 1 ;
			$name =~ s/^\s+// ;
			$seq = '' ;
		} else {
			$seq .= $_ ;
		}
	}
	$ret{$name} = $seq if $seq ne '' ;
	close FILE ;
	return \%ret ;
}

1 ;