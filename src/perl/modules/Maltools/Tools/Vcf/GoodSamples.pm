package Maltools::Tools::Vcf::GoodSamples;
use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  vc_qthr => { type => 'num', default => '30' },
  vt_qthr => { type => 'num', default => '10' },
  vc_qthr_hetero => { type => 'num', default => '30' },
  vt_dthr => { type => 'num', default => 5 },
  svc_cvg_thr => { type => 'num', default => 0.5 },
};

our $OUTPUTS = {
  out => {type => 'file', mandatory => 1},
};

sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
  return "perl";
}


sub command_template {
   return TTS('{$T->interpreter} {$S}');
}


sub script_template {
  return Text::Template->new(TYPE=> 'FILEHANDLE', SOURCE => \*DATA);
}

__DATA__
#!/usr/bin/env perl

use strict;
use warnings;

use IO::File;
use Vcf;

my $input = "{ $J->input("in") || '-' }";
my $output = "{ $J->output("out") || '-' }";
my $vc_qthr = { $J->input("vc_qthr") || 'undef' };
my $vt_qthr = { $J->input("vt_qthr") || 'undef' };
my $vc_qthr_hetero = { $J->input("vc_qthr_hetero") || 'undef' };
my $vt_dthr = { $J->input("vt_dthr") || 'undef' };
my $svc_cvg_thr = { $J->input("svc_cvg_thr") || '0.5' };

$svc_cvg_thr or die "no valid credible variant coverage threshold provided";

my $input_fh = $input eq "-" ? \*STDIN : IO::File->new($input,'r') or die "could not open input file '$input' for reading";
my $output_fh = $output eq "-" ? \*STDOUT : IO::File->new($output,'w') or die "could not open output file '$output' for writing";

my $vcf = Vcf->new(fh => $input_fh);

$vcf->parse_header;

my @sample_names = $vcf->get_samples;

my %sample_data = map \{ $_ => \{ name => $_, vc_count => 0, vc_count_hetero => 0, tvc_count => 0, tvc_count_hetero => 0 \} \} @sample_names;
my $vc_count = 0;
my $v_count = 0;
while (my $variant = $vcf->next_data_hash) \{
  $v_count++;
  my $qual = $variant->\{QUAL\};
  my $filters = join(";",@\{$variant->\{FILTER\}\});
  print STDERR join("\t",$variant->\{CHROM\},$variant->\{POS\},$filters),"\n";
  next if $vc_qthr && $qual < $vc_qthr;
  next unless $filters eq "" || $filters eq "PASS" || $filters eq ".";
  $vc_count++;
  process_variant($variant,\%sample_data);
\}

if ($v_count == 0) \{
  die "there is no variants to analyze in the input";
\}

if ($vc_count == 0) \{
  die "there is no credible variants present in the input";
\}

print $output_fh "##Processed $v_count variants of which $vc_count are considered credible (" . sprintf("%.2f",100 * $vc_count / $v_count) . "%)\n"; 
print $output_fh "##vc_qthr = { $J->input("vc_qthr") || 'None' }, vt_qthr = { $J->input("vt_qthr") || 'None' }, vc_qthr_hetero = { $J->input("vc_qthr_hetero") || 'None' }" . 
                 ", vt_dthr = { $J->input("vt_dthr") || 'None' }, svc_cvg_thr = { $J->input("svc_cvg_thr") || '0.5' }\n";
print $output_fh "##Sample summary:\n";
print $output_fh "#SAMPLE\tCLASS\tTYPABLE_COVERAGE\tTYPABLE\tTYPABLE_HETERO\tCREDIBLE_COUNT\tCREDIBLE_HETERO\tTYPABLE_HETEROZYGOSITY\tCREDIBLE_HETEROZYGOSITY\n";
foreach my $sname (@sample_names) \{
  my $sdata = $sample_data\{$sname\};
  $sdata->\{vc_cvg\} = $sdata->\{vc_count\} / $vc_count;
  $sdata->\{class\} = $sdata->\{vc_cvg\} < $svc_cvg_thr ? "POOR" : "GOOD";
  print $output_fh join("\t",$sname,$sdata->\{class\},$sdata->\{vc_cvg\},$sdata->\{vc_count\},$sdata->\{vc_count_hetero\},$sdata->\{tvc_count\},
          $sdata->\{tvc_count_hetero\},
          $sdata->\{tvc_count\} ? $sdata->\{tvc_count_hetero\}/$sdata->\{tvc_count\} : "NA",
          $sdata->\{vc_count\} ? $sdata->\{vc_count_hetero\}/$sdata->\{vc_count\} : "NA" ),"\n";
\}

$output_fh->close() unless $output eq "-";


sub process_variant \{
  my ($variant,$sd_hash) = @_;
  foreach my $sname (keys %$sd_hash) \{
    my $gtype = $variant->\{gtypes\}->\{$sname\};
    my $gq = $gtype->\{GQ\} || '0';
    my $dp = $gtype->\{DP\} || '0';
    my $gt = $gtype->\{GT\} || './.';
    next if $gt =~ /\..?\./; # disregard deletions.
    #print STDERR "$sname $gq $dp $gt\n";
    my ($a1,$a2) = split(/\/|\|\\/,$gt);
    my $hetero = $a1 ne '.' && $a2 ne '.' && $a1 != $a2;
    $sd_hash->\{$sname\}->\{tvc_count\}++;
    $sd_hash->\{$sname\}->\{tvc_count_hetero\}++ if $hetero;
    next if $vt_qthr && $gq < $vt_qthr;
    next if $vt_dthr && $dp < $vt_dthr;
    next if $vc_qthr_hetero && $hetero && $gq < $vc_qthr_hetero;    
    $sd_hash->\{$sname\}->\{vc_count_hetero\}++ if $hetero;
    $sd_hash->\{$sname\}->\{vc_count\}++;
  \}
\}
