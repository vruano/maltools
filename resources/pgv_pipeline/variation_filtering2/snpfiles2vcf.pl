#!/usr/bin/env perl

# NOTE : Runtime ~4h for 425 samples/1M SNPs, not counting the sorted tmp file generation (another 1-2h)

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 0 ; # Set to 1 to clean up temporary files
my $show_gt = 1 ;

my ( $manifest_file , $sample_range ) ;
my $batch_size ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'batchsize=i' => \$batch_size ,
				'samplerange=s' => \$sample_range # Use samplerange to process only a subset of samples in the manifest, then join the VCF files later. Untested. 1-based.
) ;

die "Need --manifest!\n" unless defined $manifest_file ;

my $is_farm = ( defined $ENV{'LSB_JOBINDEX'} and defined $batch_size ) ? 1 : 0 ;
my ( $sample_range_from , $sample_range_to ) ;
if ( $is_farm ) {
	$sample_range_from = ( $ENV{'LSB_JOBINDEX'} - 1 ) * $batch_size + 1 ; # 1-based
	$sample_range_to = ( $ENV{'LSB_JOBINDEX'} ) * $batch_size - 1 + 1 ; # 1-based
} elsif ( defined $sample_range ) {
	$sample_range =~ m/^(\d+)-(\d+)$/ ;
	( $sample_range_from , $sample_range_to ) = ( $1 , $2 ) ;
}

# Prep
read_manifest ( $manifest_file ) ;
my $snpfile_dir = $manifest->{'paths'}->{'snpdata_path'} || die "No snpdata_path in manifest file\n" ;
my $tmpdir ;
$tmpdir = "$outdir/tmp.snpfiles2vcf" ;
$tmpdir .= ".$sample_range_from-$sample_range_to" if defined $sample_range_from ;
`rm -rf $tmpdir` if $production and -d $tmpdir ; # Activate for production; ensures clean temporary data
mkdir $tmpdir unless -d $tmpdir ;
my $sort = 'sort -k1,1 -k2,2n' ;

# Read SNP data
my %snps ;
print STDERR "Reading SNP list\n" unless $production or $is_farm ;
open SNPLIST , get_manifest_path ( 'candidate_snps' ) ;
while ( <SNPLIST> ) {
    next unless $_ =~ m/^(\S+)\s(\d+)\s(.)\s(.)$/ ;
    $snps{$1}->{$2} = [ $3 , $4 , 99 , 0 ] ; # ref, alt, uniqueness, coding
}
close SNPLIST ;

# Read uniqueness data
open UNIQUENESS , get_manifest_path ( 'uniqueness_scores' ) ;
while ( <UNIQUENESS> ) {
	next unless $_ =~ m/^(\S+)\s(\d+)\s(\d+)$/ ;
	next unless defined $snps{$1}->{$2} ;
	$snps{$1}->{$2}->[2] = $3 ;
}
close UNIQUENESS ;

# Read SNP properties
open SNPPROPS , get_manifest_path ( 'snp_properties' ) ;
while ( <SNPPROPS> ) {
	next unless $_ =~ m/^(\S+)\s(\d+)\s\S+\s\S+\s(\S+)\s/ ;
	next unless defined $snps{$1}->{$2} ;
	my ( $chr , $pos , $type ) = ( $1 , $2 , $3 ) ;
	$snps{$chr}->{$pos}->[3] = 1 if $type =~ m/synonymous$/ ;
}
close SNPPROPS ;


# Create VCF object
my $vcf = Vcf->new ( version => '4.0' ) ; # file => $vcf_file , 
$vcf->recalc_ac_an ( 0 ) ;

# Sort SNP files into temporary directory
my %filenames ;
my $sample_count = 0 ;
foreach my $sample ( @{$manifest->{'samples'}} ) {
	$sample_count++ ;
	next if defined $sample_range_from and $sample_count < $sample_range_from ;
	last if defined $sample_range_to and $sample_count > $sample_range_to ;
    my $infile = "$snpfile_dir/$sample.snps" ;
    unless ( -e $infile ) {
        print STDERR "Skipping missing SNP data file for $sample : $infile\n" ;
        next ;
    }
    if ( 0 == -s $infile ) {
        print STDERR "Skipping empty SNP data file for $sample : $infile\n" ;
        next ;
    }
	my $outfile = "$tmpdir/$sample.tmp" ;
	$filenames{$sample} = $outfile ;
	run ( 'gawk \'{ print $2 "\t" $3 "\t" "' . $sample . '" "\t" $1 "\t" $4 "\t" $5 }\' ' . $infile . " | $sort > " . $outfile ) unless -e $outfile ;
	$vcf->add_columns ( $sample ) ;
}

# The dummy sample ensures that all SNPs from the candidate list are listed, even if they do not have a single call; it will not show up in the VCF file
my $dummysample = "osuhrvsi7hsvio7shils4gvkdrdv" ; # Some hopefully unique name
my $outfile = "$tmpdir/$dummysample.tmp" ;
$filenames{$dummysample} = $outfile ;
unless ( -e $outfile ) {
	print STDERR "Creating dummy sample\n" unless $production or $is_farm ;
    open OUT , '| gawk \'{ print $2 "\t" $3 "\t" "' . $dummysample . '" "\t" $1 "\t" $4 "\t" $5 }\' ' . " | $sort > $outfile " ;
    foreach my $chr ( keys %snps ) {
        foreach my $pos ( keys %{$snps{$chr}} ) {
            print OUT "1\t$chr\t$pos\t" . $snps{$chr}->{$pos}->[0] . "\t" . $snps{$chr}->{$pos}->[0] . "\n" ;
        }
    }
    close OUT ;
}


# Header
$vcf->add_header_line ( {key=>'INFO', ID=>'AD',Number=>'-1',Type=>'Integer',Description=>'Allele count in genotypes'} ) ;
$vcf->add_header_line ( {key=>'INFO', ID=>'CODING',Number=>'0',Type=>'Flag',Description=>'Mark variants at protein coding positions'} ) ;
$vcf->add_header_line ( {key=>'INFO', ID=>'UQ',Number=>'1',Type=>'Integer',Description=>'Uniqueness score for this variant'} ) ;
$vcf->add_header_line ( {key=>'INFO', ID=>'NS',Number=>'1',Type=>'Integer',Description=>'Number of Samples With Data'} ) ;
$vcf->add_header_line ( {key=>'INFO', ID=>'DP',Number=>'1',Type=>'Integer',Description=>'Total Depth'} ) ;
#$vcf->add_header_line ( {key=>'INFO', ID=>'AN',Number=>'1',Type=>'Integer',Description=>'Total number of alleles in called genotypes'} ) ;

$vcf->add_header_line ( {key=>'FORMAT', ID=>'GT',Number=>1,Type=>'String',Description=>'Genotype'} ) ;
$vcf->add_header_line ( {key=>'FORMAT', ID=>'DP',Number=>1,Type=>'Integer',Description=>'Read depth'} ) ;
$vcf->add_header_line ( {key=>'FORMAT', ID=>'AD',Number=>'.',Type=>'Integer',Description=>'Allelic depths for the ref and alt alleles in the order listed, base-quality filtered'} ) ;

my $outfh ;
if ( $is_farm ) {
	open $outfh , ">$outdir/batch.".$ENV{'LSB_JOBINDEX'}.".vcf" ;
	print $outfh $vcf->format_header();
} else {
	print $vcf->format_header();
}

# Calls
my %calls ;
$calls{'chr'} = '' ;
$calls{'pos'} = '' ;
open FILE , "$sort -m $tmpdir/*.tmp |" ;
while ( <FILE> ) {
	next unless $_ =~ m/^(\S+)\s(\d+)\s(\S+)\s(\d+)\s(\S+)\s(\S+)$/ ;
	my ( $chr , $pos , $sample , $count , $observed , $reference ) = ( $1 , $2 , $3 , $4 , $5 , $6 ) ;
	write_line ( $chr , $pos ) if $chr ne $calls{'chr'} or $pos != $calls{'pos'} ;
	next if $sample eq $dummysample ;
	$calls{'data'}->{$sample}->{$observed} = $count ;
}
close FILE ;
write_line ( '' , '' ) ; # Flush last SNP


# Cleanup
`rm -rf $tmpdir` if $production ;

0 ; # This is the end, my friend, the end

sub write_line {
	my ( $chr , $pos ) = @_ ;
	
	if ( $calls{'chr'} ne '' and defined $snps{$calls{'chr'}} and defined $snps{$calls{'chr'}}->{$calls{'pos'}} ) {
		my (@samples) = $vcf->get_samples();
		my $snp = $snps{$calls{'chr'}}->{$calls{'pos'}} ;
		
		my %x ;
#		$x{'ID'} = $calls{'chr'} . ':' . $calls{'pos'} ;
		$x{'CHROM'} = $calls{'chr'} ;
		$x{'POS'} = $calls{'pos'} ;
		$x{'REF'} = $snp->[0] ;
		
		my $alt = $iupac_reverse{$snp->[1]} ;
		$alt =~ s/($snp->[0])//g ; # Remove reference
		my @bases = sort ( split ( '' , $alt ) ) ;
		my @alt = @bases ;
		unshift @bases , $snp->[0] ;
		
		$x{'FILTER'} = ['.'] ;
		$x{'QUAL'} = -1 ;
		
		$x{'INFO'} = { DP => 0 , NS => 0 , UQ => $snp->[2] } ;
		$x{'INFO'}->{'CODING'} = 1 if $snp->[3] ; # How to set this without a value? undef doesn't seem to work, neither does ''. See ~40 lines below for hack-around...
		
		my @ac ;
		push @ac , 0 foreach @bases ;
		
		$vcf->add_format_field ( \%x , 'GT' ) if $show_gt ;
		$vcf->add_format_field ( \%x , 'DP' ) ;
		$vcf->add_format_field ( \%x , 'AD' ) ;
		
		$x{'gtypes'} = {} ;
		foreach my $sample ( @samples ) {
			my $d = $calls{'data'}->{$sample} ;
			$x{'gtypes'}->{$sample}->{'GT'} = '.' if $show_gt ;
			
			$x{'gtypes'}->{$sample}->{'DP'} = 0 ;
			$x{'gtypes'}->{$sample}->{'AD'} = [] ;

			if ( defined $d ) {
				$x{'gtypes'}->{$sample}->{'DP'} += $_ foreach values %{$d} ;
				push @{$x{'gtypes'}->{$sample}->{'AD'}} , ($d->{$_}||0) foreach @bases ;
				
				foreach my $n ( 0 .. $#bases ) {
					$ac[$n] += ($d->{$bases[$n]}||0) ;
#					push @{$x{'gtypes'}->{$sample}->{'GT'}} , ( defined $d->{$bases[$n]} ? $bases[$n] : '.' ) ;
				}
				$x{'INFO'}->{'NS'}++ ;
			} else {
				push @{$x{'gtypes'}->{$sample}->{'AD'}} , 0 foreach @bases ;
#				push @{$x{'gtypes'}->{$sample}->{'GT'}} , '.' foreach @bases ;
			}
			
#			$x{'gtypes'}->{$sample}->{'GT'} = join '/' , @{$x{'gtypes'}->{$sample}->{'GT'}} ;
			$x{'gtypes'}->{$sample}->{'AD'} = join ',' , @{$x{'gtypes'}->{$sample}->{'AD'}} ;
			
			$x{'INFO'}->{'DP'} += $x{'gtypes'}->{$sample}->{'DP'} ;
		}
		
		
		$x{'ALT'} = \@alt ; # Manual allele order
		$x{'INFO'}->{'AD'} = join "," , @ac ;
		my $s = $vcf->format_line(\%x);
		$s =~ s/\bCODING=1\b/CODING/ ; # ARGH MY EYES!!!

		if ( $is_farm ) {
			print $outfh $s ;
		} else {
			print $s ;
		}
	}	
	
	my %blank ;
	$calls{'chr'} = $chr ;
	$calls{'pos'} = $pos ;
	$calls{'data'} = \%blank ;
}

sub run {
    my ( $cmd ) = @_ ;
    print STDERR "$cmd\n" unless $production or $is_farm ;
    `$cmd` ;
}
