#!/usr/bin/env perl

use strict ;
use warnings ;

BEGIN { unshift @INC , "../project_snapshot" }
use gff2json ;
use Variation ;

my ( $manifest_file , $coding_regions_only ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'coding_only' => \$coding_regions_only
) ;

die "Needs --manifest\n" unless defined $manifest_file  ;
read_manifest ( $manifest_file ) ;

my $f_reference = get_manifest_path('reference') ;
my $f_annotation = get_manifest_path('annotation');
my $f_snplist = get_manifest_path('candidate_snps');


my $fasta = read_fasta ( $f_reference ) ;
my $gff = get_gff_object_from_file ( $f_annotation ) ;
my $snps = read_snps ( $f_snplist ) ;

# Paranoia
compare_chromosome_names ( $fasta , $gff , "FASTA" , "GFF" ) ;
compare_chromosome_names ( $fasta , $snps , "FASTA" , "SNPS" ) ;
compare_chromosome_names ( $gff , $snps , "GFF" , "SNPS" ) ;

foreach my $chr ( sort keys %{$gff} ) {
    next unless defined $fasta->{$chr} ;
    next unless defined $snps->{$chr} ;
    my %genes_with_mrna ;
    foreach my $mrna ( @{$gff->{$chr}->{'cds'}} ) {
        my $rna = $gff->{$chr}->{'data'}->{$mrna}->[0]->{'parent'} ;
        my $gene = $gff->{$chr}->{'data'}->{$rna}->[0]->{'parent'} ;
        push @{$genes_with_mrna{$gene}} , $_ foreach ( @{$gff->{$chr}->{'data'}->{$mrna}} ) ;
    }
    
    foreach my $gene ( @{$gff->{$chr}->{'gene'}} ) {
        print STDERR "$gene\n" ;
        next unless defined $genes_with_mrna{$gene} ;
        @{$genes_with_mrna{$gene}} = sort { $a->{'from'} <=> $b->{'from'} } @{$genes_with_mrna{$gene}} ;
        my $desc = $gff->{$chr}->{'data'}->{$gene}->[0]->{'description'} || '-' ;
        next if $desc =~ m/\bpseudogene\b/i ;
        
        my $dna = get_coding_dna ( $gene , $chr , $genes_with_mrna{$gene} ) || '' ;

        unless ( $dna =~ m/^ATG/ ) { # About 35 of those in V5.5
            print STDERR "Gene $gene ($desc) on $chr does not start with ATG : " . substr ( $dna.'                    ' , 1, 20 ) . "\n" ;
            annotate_with_blank_gene ( $chr , $gene , $genes_with_mrna{$gene} ) ;
            next ;
        }
#        unless ( ( length $dna ) % 3 == 0 ) {
 #           print STDERR "Gene $gene on $chr has length non-dividable by 3 : " . length ( $dna ) . "\n" ;
  #          next ;
   #     }

        my $orig_aa = dna2aa ( $dna ) ;
        if ( $orig_aa eq '' ) {
            print STDERR "No amino acid sequence : $gene\t$desc\n" ;
            annotate_with_blank_gene ( $chr , $gene , $genes_with_mrna{$gene} ) ;
            next ;
        }
        
 #       print "$chr\t$gene\n$dna\n$original_aa\n\n" ;
#        print "OK : $chr\t$gene\n$desc\n" ;
        
        foreach my $mrna ( @{$genes_with_mrna{$gene}} ) {
        
        	if ( $coding_regions_only ) {
        		print "$chr\t" . $mrna->{'from'} . "\t" . $mrna->{'to'} . "\n" ;
        		next ;
        	}
        	
        	
            foreach my $pos ( $mrna->{'from'} .. $mrna->{'to'} ) {
                next unless defined $snps->{$chr}->{$pos} ;
                
                my $orig_base = substr ( $fasta->{$chr} , $pos-1 , 1 ) ;
                my @alt = split '' , $iupac_reverse{$snps->{$chr}->{$pos}->[1]} ;
                foreach my $alt ( @alt ) {
                    next if $alt eq $orig_base ;
                    substr ( $fasta->{$chr} , $pos-1 , 1 ) = $alt ;
                    my $new_dna = get_coding_dna ( $gene , $chr , $genes_with_mrna{$gene} ) || '' ;
                    my $new_aa = dna2aa ( $new_dna ) ;
                    my $syn ;
                    my $subst = '-' ;
                    if ( $orig_aa eq $new_aa ) { # Synonymous
                        $syn = 'synonymous' ;
                    } else { # Non-synonymous
                        $syn = 'nonsynonymous' ;
                        foreach ( 0 .. length ( $orig_aa ) - 1 ) {
                            next if substr($orig_aa,$_,1) eq substr($new_aa,$_,1) ;
                            $subst = substr($orig_aa,$_,1) . ($_+1) . substr($new_aa,$_,1) ;
                            last ;
                        }
                    }
                    my $s = get_row ( $chr , $pos , $alt , $syn , $subst , $gene ) ;
#                    print STDERR "$s\n" ;
                    push @{$snps->{$chr}->{$pos}->[2]} , $s ;
                }
                substr ( $fasta->{$chr} , $pos-1 , 1 ) = $orig_base ;
                
            }
        }
    }
    
    next if $coding_regions_only ;
    
    my @pos = sort { $a <=> $b } keys %{$snps->{$chr}} ;

    # annotate SNPs in genes that have no annotation
    foreach my $gene ( @{$gff->{$chr}->{'gene'}} ) {
        foreach my $pos ( @pos ) {
            next if $pos < $gff->{$chr}->{'data'}->{$gene}->[0]->{'from'} ;
            last if $pos > $gff->{$chr}->{'data'}->{$gene}->[0]->{'to'} ;
#        foreach my $pos ( $gff->{$chr}->{'data'}->{$gene}->[0]->{'from'} .. $gff->{$chr}->{'data'}->{$gene}->[0]->{'to'} ) {
#            next unless defined $snps->{$chr}->{$pos} ;
            next if defined $snps->{$chr}->{$pos}->[2] ; # Already has annotation from exons

            my $orig_base = substr ( $fasta->{$chr} , $pos-1 , 1 ) ;
            my @alt = split '' , $iupac_reverse{$snps->{$chr}->{$pos}->[1]} ;
            foreach my $alt ( @alt ) {
                next if $alt eq $orig_base ;
                push @{$snps->{$chr}->{$pos}->[2]} , get_row ( $chr , $pos , $alt , 'intron' , '-' , $gene ) ;
            }
        }
    }
    
    # annotate all other SNPs
    foreach my $pos ( keys %{$snps->{$chr}} ) {
        next if defined $snps->{$chr}->{$pos}->[2] ; # Already has annotation from exons

        my $orig_base = $snps->{$chr}->{$pos}->[0] ;
        my @alt = split '' , $iupac_reverse{$snps->{$chr}->{$pos}->[1]} ;
        foreach my $alt ( @alt ) {
            next if $alt eq $orig_base ;
            push @{$snps->{$chr}->{$pos}->[2]} , get_row ( $chr , $pos , $alt , 'intergenic' , '-' , undef ) ;
        }
    }
}

exit if $coding_regions_only ;

print "Chromosome\tPosition\tReference\tAlternate\tType\tAA_change\tGene\tGene_alt\tGene_desc\n" ;
foreach my $chr ( sort keys %{$snps} ) {
    foreach my $pos ( sort { $a <=> $b } keys %{$snps->{$chr}} ) {
        foreach ( @{$snps->{$chr}->{$pos}->[2]} ) {
            print "$_\n" ;
        }
    }
}

0 ;

sub annotate_with_blank_gene {
    my ( $chr , $gene , $mrnas ) = @_ ;
    foreach my $mrna ( @{$mrnas} ) {
        foreach my $pos ( $mrna->{'from'} .. $mrna->{'to'} ) {
            next unless defined $snps->{$chr}->{$pos} ;

            my $orig_base = substr ( $fasta->{$chr} , $pos-1 , 1 ) ;
            my @alt = split '' , $iupac_reverse{$snps->{$chr}->{$pos}->[1]} ;
            foreach my $alt ( @alt ) {
                next if $alt eq $orig_base ;
                push @{$snps->{$chr}->{$pos}->[2]} , get_row ( $chr , $pos , $alt , 'unknown' , '-' , $gene ) ;
            }

        }
    }
}

sub get_row {
    my ( $chr , $pos , $alt , $syn , $subst , $gene ) = @_ ;
    my $orig_base = $snps->{$chr}->{$pos}->[0] ;
    return "$chr\t$pos\t$orig_base\t$alt\t$syn\t$subst\t" . get_gene_desc ( $chr , $gene ) ;
}

sub get_gene_desc {
    my ( $chr , $gene ) = @_ ;
    return "-\t-\t-" unless defined $gene ;
    my $desc = $gff->{$chr}->{'data'}->{$gene}->[0]->{'description'} || '-' ;
    my $alt = '-' ; # FIXME
    return "$gene\t$alt\t$desc" ;
}

sub dna2aa {
    my ( $dna ) = @_ ;
    my $aa ;
    while ( $dna =~ m/(...)/g ) {
        $aa .= $codon2aa{$1} || 'X' ;
#        last if 'X' eq ( $codon2aa{$1} || 'X' ) ;
    }
    return '' unless $aa =~ m/^([^X]*)X(.*)$/ ; # Stop codon
#    return '' unless '' eq $2 ;
    $a = $1 || '' ;
    return $aa ;
}

sub get_coding_dna {
    my ( $gene , $chr , $mrnas ) = @_ ;
    my @p2 ;
    my $strand ;
    foreach my $mrna ( @{$mrnas} ) {
        $strand = $mrna->{'strand'} ;
        push @p2 , $_ foreach ( $mrna->{'from'} .. $mrna->{'to'} ) ;
#        print "mRNA : " . $mrna->{'from'} . "-" . $mrna->{'to'} . "\n" ;
    }
    unless ( defined $strand ) { # No strand defined
        print STDERR "No strand defined for $gene\n" ;
        return '' ;
    }
    
    my %p ;
    $p{$_} = 1 foreach @p2 ;
    my @p = sort { $a <=> $b } keys %p ;
    if ( scalar @p != scalar @p2 ) {
        print STDERR "Overlapping exons in $gene\n" ;
        return '' ;
    }
    @p = reverse @p if $strand eq '-' ;
    
    my $ret = '' ;
    $ret .= substr ( $fasta->{$chr} , $_-1 , 1 ) foreach @p ;
    $ret =~ tr/ACGT/TGCA/ if $strand eq '-' ;
    return $ret ;
}

sub read_snps {
    my ( $file ) = @_ ;
    my %ret ;
    open FILE , $file ;
    while ( <FILE> ) {
        next unless $_ =~ m/^(\S+)\s(\d+)\s(.)\s(.)$/ ;
        $ret{$1}->{$2} = [ $3 , $4 ] ;
        die "Base mismatch between SNP list and reference at $1:$2\n" unless $3 eq substr $fasta->{$1} , $2 - 1 , 1 ;
    }
    close FILE ;
    return \%ret ;
}

sub compare_chromosome_names {
    my ( $h1 , $h2 , $hn1 , $hn2 , $not_again ) = @_ ;
    foreach ( sort keys %{$h1} ) {
        next if defined $h2->{$_} ;
        print STDERR "ATTENTION : Chromosome $_ defined in $hn1, but not in $hn2\n" ;
    }
    compare_chromosome_names ( $h2 , $h1 , $hn2 , $hn1 , 1 ) unless $not_again ;
}

