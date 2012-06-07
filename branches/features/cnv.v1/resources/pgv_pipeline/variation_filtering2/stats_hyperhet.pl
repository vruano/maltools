#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;
use Vcf ;

my $production = 1 ;

my ( $manifest_file , $vcf_file , $group , $update_plot_only ) ;
my $amin = 2 ;
my $cmin = 5 ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'vcf=s' => \$vcf_file , 
				'group=s' => \$group , 
				'amin=i' => \$amin ,
				'cmin=i' => \$cmin ,
				'onlyplot' => \$update_plot_only
) ;

my $usage = "USAGE : ./stats_hyperhet.pl --manifest MANIFEST_FILE --vcf FILE --group GROUP_NAME [--amin INT|2] [--cmin INT|5] [--onlyplot]\n" ;
die $usage unless defined $manifest_file ;
die $usage unless defined $vcf_file ;
die $usage unless defined $group ;

read_manifest ( $manifest_file ) ;

if ( $group eq 'all' ) {
	my @groups = sort keys %{$manifest->{'groups'}} ;
	
	my $cmd = '' ;
	foreach ( @groups ) {
		my $file = $vcf_file ;
		$file = '-' unless $cmd eq '' ; # Pipe on all but first
		$cmd .= "./stats_hyperhet.pl --manifest $manifest_file --vcf $file --group $_ --amin $amin --cmin $cmin |" ;
	}

	open FILE , $cmd ;
	while ( <FILE> ) {
		print $_ ;
	}
	close FILE ;
	
	exit 0 ;
}

die "Group $group not defined in manifest $manifest_file\n" unless defined $manifest->{'groups'}->{$group} ;
die "Sample list for $group not defined in manifest $manifest_file\n" unless defined $manifest->{'groups'}->{$group}->{'samples'} ;

my $vcf = Vcf->new ( file => $vcf_file , ignore_missing_GT => 1 ) ;
$vcf->recalc_ac_an ( 0 ) ;
$vcf->parse_header();

my (@current_samples) = $vcf->get_samples();

my %group_samples ;
$group_samples{$_} = 1 foreach ( @{$manifest->{'groups'}->{$group}->{'samples'}} ) ;

my @group_sample_ids ; # Will hold array indices for group samples
foreach my $sid ( 0 .. $#current_samples ) {
	my $sample = $current_samples[$sid] ;
	push @group_sample_ids , $sid if defined $group_samples{$sample} ;
}

my $out_base = $manifest->{'paths'}->{'out_path'} . '/hyperheterozygous' ;
mkdir $out_base unless -d $out_base ;
$out_base .= "/hyperhet.$group" ;
my $out_pls  = "$out_base.pls" ;

my @freq ;
push @freq , "Chr\tPos\tRef\tAlleles\tA\tB\trA\trB\tfA\tfB\tfMAF\tBin\tph\tt\th" ;

print $vcf->format_header() if $production ;

my $ad ;
my $dp ;
my @alleles ;

# Parse through data
while (my $x=$vcf->next_data_array()) {

	if ( $x->[6] eq 'PASS' ) {
	
		my @o ;
		foreach ( @group_sample_ids ) {
			my @m = split ( ':' , $x->[9+$_] ) ;
			push @o , \@m ; 
		}
		
		$ad = get_format_id ( $x , 'AD' ) ;
		$dp = get_format_id ( $x , 'DP' ) ;
		@alleles = ( $x->[3] , split ( ',' , $x->[4] ) ) ;
		
		push @freq , generate_frequency ( $x , \@o ) ;
	}
	
	print join ( "\t" , @{$x} ) . "\n" if $production ;
}

$vcf->close() ;

calculate_pls ( $out_pls ) unless $update_plot_only and -e $out_pls ;
generate_R_plot ( $out_pls , $out_base , $manifest->{'groups'}->{$group}->{'hyperhet_cutoff'} ) ;

0 ;
# END OF MAIN


sub get_format_id {
	my ( $x , $key ) = @_ ;
	my @format = split ':' , $x->[8] ;
	my %format ;
	$format{$format[$_]} = $_ foreach ( 0 .. $#format ) ;
	return $format{$key} ;
}

sub generate_frequency {
	my ( $x , $o ) = @_ ;

	my $chr = $x->[0] ;
	my $pos = $x->[1] ;
	my $ref = $x->[3] ;
	
	my $count_t = 0 ;
	my $count_h = 0 ;
	my @arr_t ;
	my @arr_h ;
	my %r ;
	$r{'CHR'} = $chr ;
	$r{'POS'} = $pos ;
	$r{'REF'} = $ref ;
	foreach my $d ( @{$o} ) {
		my $sum = 0 ;
		my $above_amin = 0 ;
		
		my @counts = split ',' , $d->[$ad] ;
		foreach my $p ( 0 .. $#counts ) {
			my $base = $alleles[$p] ;
			my $count = $counts[$p] ;
			next if $count == 0 ;
			$r{'A'} = $base unless defined $r{'A'} ;
			$r{'B'} = $base if $r{'A'} ne $base and not defined $r{'B'} ;
			if ( $r{'A'} eq $base ) {
				$r{'rA'} += $count ;
			} elsif ( $r{'B'} eq $base ) {
				$r{'rB'} += $count ;
			} else {
				die "Third base!\n" ;
			}
			$sum += $count ;
			$above_amin++ if $count >= $amin ;
		}
		
		if ( $sum >= $cmin ) {
			$count_t++ ;
			push @arr_t , 1 ;
			$count_h++ if $above_amin == 2 ;
			push @arr_h , ( $above_amin == 2 ) ? 1 : 0 ;
		} else {
			push @arr_t , 0 ;
			push @arr_h , 0 ;
		}
	}
	$r{'ALLELES'} = 0 ;
	$r{'ALLELES'}++ if defined $r{'A'} ;
	$r{'ALLELES'}++ if defined $r{'B'} ;
	if ( defined $r{'B'} ) {
		die "Hidden third base at $chr:$pos!\n" if $ref ne $r{'A'} and $ref ne $r{'B'} ; # Paranoia
		$r{'fA'} = $r{'rA'} / ( $r{'rA'} + $r{'rB'} ) ;
		$r{'fB'} = $r{'rB'} / ( $r{'rA'} + $r{'rB'} ) ;
		$r{'fMAF'} = $r{'fA'} < $r{'fB'} ? $r{'fA'} : $r{'fB'} ;
		$r{'BIN'} = 0 ;
		$r{'BIN'}++ while ( $r{'fMAF'} > ($r{'BIN'}+1) * 0.05 ) ;
		$r{'ph'} = ( $count_t > 0 ) ? ( $count_h / $count_t ) : 0 ;
		$r{'T'} = join '' , @arr_t ;
		$r{'H'} = join '' , @arr_h ;
	} else {
		$r{'A'} = '' unless defined $r{'A'} ;
		$r{'B'} = '' unless defined $r{'B'} ;
	}
	
	my $ret = '' ;
	$ret .= sprintf "%s\t%d\t%s\t%d\t%s\t%s" , $chr , $pos , $ref , $r{'ALLELES'} , $r{'A'} , $r{'B'} ;
	$ret .= sprintf "\t%d\t%d\t%4.4f\t%4.4f\t%4.4f\t%d\t%4.4f\t%s\t%s" , $r{'rA'} , $r{'rB'} , $r{'fA'} , $r{'fB'} , $r{'fMAF'} , $r{'BIN'} , $r{'ph'} , $r{'T'} , $r{'H'} if $r{'ALLELES'} == 2 ;
	return $ret ;
}

sub calculate_pls {
	my ( $out_file ) = @_ ;
	
	# Read bin data
	my @phb ;
	foreach ( @freq ) {
		my ( $chr , $pos , $ref , $alleles , $a , $b , $ra , $rb , $fa , $fb , $fmaf , $bin , $ph ) = split "\t" , $_ ;
		next if lc $chr eq 'chr' ; # Header line
		next unless $alleles == 2 ; # Only biallelic SNPs
		push @{$phb[$bin]} , $ph ;
	}
	
	# Calculate bin probabilities
	foreach my $bin ( 0 .. $#phb ) {
		$phb[$bin] = [] unless defined $phb[$bin] ; # Paranoia
		my $sum = 0 ;
		$sum += $_ foreach ( @{$phb[$bin]} ) ;
		my $n = scalar @{$phb[$bin]} ;
		$phb[$bin] = 0 == $n ? 0 : $sum / $n ;
	}
	
	# Read freqs and write pls
	open OUT , ">$out_file" ;
	foreach ( @freq ) {
		my $line = $_ ;
		my ( $chr , $pos , $ref , $alleles , $a , $b , $ra , $rb , $fa , $fb , $fmaf , $bin , $ph , $t , $h ) = split "\t" , $line ;
		if ( lc $chr eq 'chr' ) { # Header line
			print OUT "$line\tPLS\n" ;
			next ;
		}
		unless ( $alleles == 2 ) { # Only biallelic SNPs
			print OUT "$line\n" ;
			next ;
		}
		
		my @t = split '' , $t ;
		my @h = split '' , $h ;
		my @po ;
		foreach my $s ( 0 .. $#t ) {
			if ( $h[$s] == 1 ) {
				push @po , $phb[$bin] ;
			} elsif ( $h[$s] == 0 and $t[$s] == 1 ) {
				push @po , 1 - $phb[$bin] ;
			} else {
				push @po , 1 ;
			}
		}
		
		my $p = 1 ;
		$p *= $_ foreach ( @po ) ;
		my $pls = - log($p) / log(10) ;
		
		print OUT "$line\t$pls\n" ;
	}
	close OUT ;
}

sub generate_R_plot {
	my ( $pls_file , $outbase , $cutoff ) = @_ ;
	
	my $png_file = "$outbase.png" ;
	my $tmp_file = "$outbase.tmp" ;

	my $cmd = "gawk '{ if ( \$4 != 1 ) print \$16 }' $pls_file > $tmp_file" ;
	`$cmd` ;
	
	open FILE , ">$tmp_file.R" ;
	print FILE "stick_to_console='yes'\n" ;
	print FILE "d <- read.table ( \"$tmp_file\",header=T,sep=\"\t\")\n" ;
	print FILE "d2 <- sort(d[,1])\n" ;
	print FILE "bitmap(file=\"$png_file\", width=1280, height=1024, type=\"png16\", units=\"px\")\n" ;
	print FILE "plot(d2)\n" ;
	print FILE "abline(h=$cutoff,col=\"red\")\n" if defined $cutoff ;
	print FILE "axis(1,tck=1,col = \"grey\", lty = \"dotted\")\n" ;
	print FILE "axis(2,tck=1,col = \"grey\", lty = \"dotted\")\n" ;
	print FILE "dev.off()\n" ;
	close FILE ;
	`R --vanilla < $tmp_file.R` ;
	
	unlink "$tmp_file" ;
	unlink "$tmp_file.R" ;
}
