package Maltools::Tools::PileupToVarLists;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use vars qw(%ENV);
use File::Basename qw(dirname);
use Text::Template;
use IO::File;
use Cwd qw(realpath);
use JSON::XS;
use File::Spec::Functions qw(catfile file_name_is_absolute);
use POSIX;

our $INPUTS = {
   "in" => { type => 'file', multiple => 1, mandatory => 1},
};


our $OUTPUTS = { 
   "snpout" => { type => 'file', mandatory => 1 },
   "indelout" => { type => 'file', mandatory => 0 },
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 2000;
}

sub calculate_cpu_time {
  return 4 * 60 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return 'perl';
}


1;

__DATA__
{  
 $snpout = $J->output('snpout') || '/dev/null';
 @inputs = @{$J->input('in')};
 '' }


use strict ;
use warnings ;

my $snpout = q({$snpout});
my @infiles = ('{join("','",@inputs}');

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
foreach ( keys %iupac ) \{
        my $v = $iupac\{$_\} ;
        $iupac_reverse\{$v\} = $_ ;
\}

while ($#files >= 0) \{
  my $file_name = shift @infiles;
  open(FILE,$file_name) or die "could not open '$file_name' read-only");
  while ( <FILE> ) \{
        $_ =~ /^([\S]+)\s+(\d+)\s+(.)\s+(.)\s*$/ ;
        next unless defined $4 ;
        $read++ ;
        
        print STDERR "REFERENCE MISMATCH!\n" if ( $ref\{$1\}->[$2] || $3 ) ne $3 ;
        
        $ref\{$1\}->[$2] = $3 ;
        $alt\{$1\}->[$2] = merge_codes ( $4 , $alt\{$1\}->[$2] || $4 ) ;
  \}
  close(FILE)
\}

open(OUT,$snpout) or die "cannot open output '$snpout' to write";
my $out = 0 ;
foreach my $chr ( sort keys %alt ) \{
        my @data = @\{$alt\{$chr\}\} ;
        foreach ( 0 .. $#data ) \{
                next unless defined $data[$_] ;
                print OUT "$chr\t$_\t" . $ref\{$chr\}->[$_] . "\t" . $alt\{$chr\}->[$_] . "\n" ;
                $out++ ;
        \}
\}
close(OUT);

print STDERR "Read $read SNPs\n" ;
print STDERR "Wrote $out SNPs\n" ;

sub merge_codes \{
        my ( $n , $o ) = @_ ;
        return $n if $n eq $o ;
        my $new_all = $iupac\{$n\} . $iupac\{$o\} ;
        my $k = '' ;
        $k .= 'A' if $new_all =~ m/A/ ;
        $k .= 'C' if $new_all =~ m/C/ ;
        $k .= 'G' if $new_all =~ m/G/ ;
        $k .= 'T' if $new_all =~ m/T/ ;
        
        return $iupac_reverse\{$k\} ;
\}

