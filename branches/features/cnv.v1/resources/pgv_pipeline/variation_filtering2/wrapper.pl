#!/usr/bin/env perl

use strict ;
use warnings ;
use Variation ;

my $gz = 0 ;
my ( $manifest_file , $last_stage , $skip_existing ) ;

GetOptions  (	'manifest=s' => \$manifest_file ,
				'stage=i' => \$last_stage ,
				'skip' => \$skip_existing
) ;

my $usage = "./wrapper.pl --manifest MANIFEST_FILE [--stage LAST_STAGE] [--skip]\n" ;
die $usage unless defined $manifest_file  ;

read_manifest ( $manifest_file ) ;

# Sanity checks
die "No samples in manifest\n" unless defined $manifest->{'samples'} and 0 < scalar @{$manifest->{'samples'}} ;
die "No sequences in manifest\n" unless defined $manifest->{'sequences'} and 0 < scalar @{$manifest->{'sequences'}} ;
die "Paths not defined in manifest\n" unless defined $manifest->{'paths'} ;

# Checking path definitions
foreach ( qw ( base_path annotation candidate_snps reference uniqueness_scores snp_properties coding_regions coverage_cutoffs snpdata_path coverage_path out_path ) ) {
	die "Path for \"$_\" not defined in manifest\n" unless defined $manifest->{'paths'}->{$_} ;
}

# Checking essential files
foreach ( qw ( annotation candidate_snps reference ) ) {
	my $file = get_manifest_path ( $_ ) ;
	die "File $file (for \"$_\") does not exist\n" unless -e $file ;
	die "File $file (for \"$_\") is empty\n" unless 0 < -s $file ;
}

#_________________________________________________________________________________________________________
# Generate derived meta-data files

# Uniqueness file
# DEPENDENCIES : reference sequence
my $uniqueness_scores = get_manifest_path('uniqueness_scores') ;
unless ( -e $uniqueness_scores and 0 < -s $uniqueness_scores ) {
	my $reference = get_manifest_path('reference') ;
	print "Generating uniqueness scores at $uniqueness_scores (run will take ~4h and <6GB for Pf)\n" ;
	`./generate_uniqueness_scores $reference > $uniqueness_scores` ;
}

# SNP properties file
# DEPENDENCIES : reference sequence, candidate SNPs, GFF annotation
my $snp_properties = get_manifest_path('snp_properties') ;
unless ( -e $snp_properties and 0 < -s $snp_properties ) {
	print "Generating snp properties at $snp_properties (run will take ~4h)\n" ;
	`./generate_snp_properties.pl --manifest $manifest_file > $snp_properties` ;
}

# Coding regions file
# DEPENDENCIES : reference sequence, candidate SNPs, GFF annotation
my $coding_regions = get_manifest_path('coding_regions') ;
unless ( -e $coding_regions and 0 < -s $coding_regions ) {
	print "Generating coding regions at $coding_regions\n" ;
	`./generate_snp_properties.pl --manifest $manifest_file --coding_only > $coding_regions` ;
}

# Generate coverage cutoffs file
# DEPENDENCIES : coding regions, coverage files
my $coverage_cutoffs = get_manifest_path('coverage_cutoffs') ;
unless ( -e $coverage_cutoffs and 0 < -s $coverage_cutoffs ) {
	print "Generating coverage_cutoffs at $coverage_cutoffs (run will take ~3h per 100 samples)\n" ;
	`./generate_coverage_cutoffs.pl --manifest $manifest_file` ;
}

#_________________________________________________________________________________________________________
# Create initial VCF file
my $out = get_manifest_path('out_path');
`cp $manifest_file $out/manifest.json` ; # Copy manifest
my $original_vcf = $out . '/' . $manifest->{'paths'}->{'original_vcf'} ;
unless ( ( -e $original_vcf and 0 < -s $original_vcf ) or ( -e "$original_vcf.gz" and 0 < -s "$original_vcf.gz" ) ) {
	print "Generating original VCF file\n" ;
	my $cmd = "./snpfiles2vcf.pl --manifest $manifest_file" . pipe2file ( $original_vcf ) ;
	`$cmd` ;
}

# Processing
my $input_vcf = $original_vcf ;
my $stages = scalar @{$manifest->{'stages'}} ;

my %paths = map { $_ => get_manifest_path($_) } keys %{$manifest->{paths}};

foreach my $stage ( 1 .. $stages ) {
	last if defined $last_stage and $stage > $last_stage ;
	print "Processing stage $stage\n" ;
	my @steps = @{$manifest->{'stages'}->[$stage-1]} ;
	
	die "Both .vcf and .vcf.gz exist, can't decide which one to use\n" if -e $input_vcf and -e "$input_vcf.gz" ;

	my $cmd = '' ;
	if ( -e "$input_vcf.gz" ) {
		$input_vcf = "$input_vcf.gz" ;
		$cmd = "bgzip -cd $input_vcf | " ;
	}
	my $names = [] ;
	foreach my $step ( @steps ) {
		if ( '' ne $cmd ) {
			$cmd .= " | " . $step->{'cmd'} . " --manifest $manifest_file --vcf -" ;
		} else {
			$cmd .= $step->{'cmd'} . " --manifest $manifest_file --vcf $input_vcf" ;
		}
		if ( defined $step->{'params'} ) {
			foreach my $k ( keys %{$step->{'params'}} ) {
				$cmd .= " --$k" ;
				my $p = $step->{'params'}->{$k} ;
				if ( $p =~ m/^\d+$/ ) {
					$cmd .= "=$p" ;
				} else {
                                    while ($p =~ s/\$\{paths\.([^\}]+)\}/$paths{$1}/g) { };
				    $cmd .= "=\"$p\"" ;
				}
			}
		}
		push @{$names} , $step->{'name'} unless $step->{'name'} =~ m/^\./ ;
	}
	my $out_vcf = "$out/" . join ( '.' , @{$names} ) . ".vcf" ;
	$cmd .= pipe2file ( $out_vcf ) ;
	
	$out_vcf = "$out_vcf.gz" if $gz ;
	unless ( $skip_existing and -e $out_vcf and 0 < -s $out_vcf ) {
		print "Running " . join ( ', ' , @{$names} ) . "\n" ;
		#print "$cmd\n" ;
		`$cmd` ;
		print "Output written to $out_vcf\n" ;
	} else {
		print "Skipping output - already exists : $out_vcf\n" ;
	}
	$input_vcf = $out_vcf ;
}

0 ;

sub pipe2file {
	my ( $file ) = @_ ;
	if ( $gz ) {
		return " | bgzip -c > $file.gz" ;
	} else {
		return " > $file" ;
	}
}
