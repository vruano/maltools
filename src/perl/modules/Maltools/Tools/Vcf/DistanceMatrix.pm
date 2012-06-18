package Maltools::Tools::Vcf::DistanceMatrix;

use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  type => { type => 'string', default => 'any' },
};

our $OUTPUTS = {
  out => {type => 'file', mandatory => 1},
  summary => {type => 'file', mandatory => 0},
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

1;

__DATA__

use strict;
use warnings;

use Vcf;
use IO::File;

my @inputs = ({ join (",", map { "\"$_\"" } $J->input("in") )});
my $output = "{ $J->output("out") }" || "-";
my $type = "any";
my $summary = "{ $J->output("summary") }";

my %samples = ();


my @sample_names;
my $sample_count;
my @zero_per_sample;

sub process_input \{
  my $input = shift;
  print STDERR "Processing $input...";

  my $fh = $input eq "-" ? \*STDIN : IO::File->new($input,"r");
  my $vcf = Vcf->new(fh => $fh);
  $vcf->recalc_ac_an(0);
  $vcf->parse_header();

  if (scalar(@sample_names) == 0) \{
    @sample_names = $vcf->get_samples();
    $sample_count = scalar(@sample_names);
    @zero_per_sample = map \{ 0 \} @sample_names;
  \}
  else \{
    my @other_samples = $vcf->get_samples();
    for (my $i = 0; $i <= $#other_samples;$i++) \{
      $other_samples[$i] eq $sample_names[$i] or die "there are different sample order across files";
    \}
  \}

  my @sample_data = map \{ \{ name => $_, total => 0, pw_total => [@zero_per_sample], pw_distance => [@zero_per_sample] \} \} @sample_names;
  while (my $x = $vcf->next_data_hash) \{
    process_snp($x,\@sample_data);
  \}
  $vcf->close();
  $fh->close();

  my $outstr = "";
  for (my $i = 0; $i < $sample_count; $i++) \{
    my $sn = $sample_names[$i];
    $sn .= " " x (10 - length($sn));
    my $sd = $sample_data[$i]; 
    my $pw_total = $sd->\{pw_total\};
    my $pw_distance = $sd->\{pw_distance\};
    my @distances = map \{ ($pw_total->[$_] == 0) ? "NA" : ($pw_distance->[$_] / $pw_total->[$_]) \} (0 .. $sample_count - 1);
    @distances = @distances[0 .. $i - 1];
    $outstr .=  $sn  . join(" ",map \{ sprintf("%1.3f", $_) \} @distances) . "\n";
  \}
  print STDERR " done\n";
  return $outstr;
\}

my $file_count = $#inputs + 1 or die "you must indicate at least one input\n";
open(STDOUT,">$output") || die "cannot open output file $output" unless $output eq "-";
my $first_file = shift(@inputs);
my $outstr = process_input($first_file);
print "  $sample_count\n";
print $outstr;
foreach my $i (@inputs) \{
 print "  $sample_count\n";
 print process_input($i); 
\}
close(STDOUT) if $output ne "-";


#print STDERR "Total per sample\n";
#print STDERR join(";", map \{ $sample_data[$_]->\{total\} \} (0 .. $sample_count -1)),"\n";

sub alleles \{
  my $gt = shift;  
  return \{ map \{ $_ => 1 \} split(/\//,$gt) \};
\}

sub covered \{
  my $al = shift;
  return scalar (grep \{ $_ ne "." \} keys %$al );
\}

sub process_snp \{
  my ($x,$sample_data) = @_;
  if ($type ne "any") \{
    my $coding = grep \{ $_ eq 'CODING' \} @\{$x->\{INFO\}\};
    return if $type ne 'coding' && $coding;
  \}
  my @sample_gtypes = map \{ $x->\{gtypes\}->\{$_\}->\{GT\} \} @sample_names;
  my @sample_alleles = map \{ alleles($_) \} @sample_gtypes;
  my @sample_covered = map \{ covered($_) \} @sample_alleles;
  for (my $i = 0; $i < $sample_count; $i++) \{
    next unless ($sample_covered[$i]);
    my $sd1 = $$sample_data[$i];
    $sd1->\{pw_total\}->[$i] = ++$sd1->\{total\};
    my $al1 = $sample_alleles[$i];
    for (my $j = $i + 1; $j < $sample_count; $j++) \{
      next unless ($sample_covered[$j]);
      my $sd2 = $$sample_data[$j];
      my $al2 = $sample_alleles[$j];
      my $intersec = scalar( grep \{ exists($al2->\{$_\}) && $_ ne "." \} keys %$al1 );
      my %union = (%$al1,%$al2);
      delete $union\{"."\};
      my $union = scalar(keys %union);
      $sd2->\{pw_total\}->[$i] = ++ $sd1->\{pw_total\}->[$j];
      $sd1->\{pw_distance\}->[$j] += ($union - $intersec) / $union; 
      $sd2->\{pw_distance\}->[$i] = $sd1->\{pw_distance\}->[$j];
    \} 
  \}
\}
