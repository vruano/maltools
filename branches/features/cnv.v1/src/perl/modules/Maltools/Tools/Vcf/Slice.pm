package Maltools::Tools::Vcf::Slice;
use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  mode => { type => 'string', mandatory => 0, default => 'last' },
  filters => { type => 'string', mandatory => 0, multiple => 1 },
  samples => { type => 'string', multiple => 1 },
  default => { type => 'string', default => 'include', multiple => 0 }, # can be 'include' or 'exclude'
};

our $OUTPUTS = {
  out => {type => 'vcf', mandatory => 1},
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
use Maltools::Vcf::Expression;
use Scalar::Util qw(blessed);

my $input = "{ $J->input("in") }";
my $output = "{ $J->output("out") }";
my $mode = "{$J->input("mode")}";
my $def = { $J->input("default") eq "include" ? 1 : 0 };
my @filters = ({join(",",map { "\"$_\"" } @{$J->input("filters")})});
my @samples = qw({join(" ",@{$J->input("samples") || []})});

open(STDIN,$input) or die "cannot open input file $input" if $input ne "-";
open(STDOUT,">$output") or die "cannot open output file $output" if $output ne "-";

my $vcf = Vcf->new(fh => \*STDIN);
$vcf->recalc_ac_an(0);
$vcf->parse_header();

@samples = $vcf->get_samples() if scalar(@samples) == 0;
@filters = map \{ process_filter($_,$vcf) \} @filters;

print $vcf->format_header(\@samples);
while (my $line = $vcf->next_line) \{
  my $x = $vcf->next_data_hash($line);
  my $fh = \{ map \{ lc($_) eq "pass" || $_ eq "." ? (PASS => 1) : ($_ => 1,NOPASS => 1) \} @\{$x->\{FILTER\}\} \};
  my $tbi = undef;
  foreach my $f (@filters) \{
    my ($name,$sign) = @$f;
    unless($name) \{
      $tbi = $sign;
    \}
    elsif (blessed($name)) \{
      $tbi = $name->evaluate($x); 
    \}
    else \{
      my $exists = exists($fh->\{$name\});
      $tbi = $sign if $exists;
    \}
    last if defined($tbi) && $name && $mode eq "first";
  \}
  $tbi = $def unless defined $tbi;
  print $vcf->format_line($x,\@samples) if $tbi; 
\}
$vcf->close();

close(STDIN) if $input ne "-";
close(STDOUT) if $output ne "-";

sub process_filter \{
  my $filter = shift;
  my $vcf = shift;
  $filter =~ s/(\+|\-)$//;
  my $sign = $1 || ($def ? '+' : '-');
  if ($filter =~ /^\s*\/(.*)\/\s*$/) \{
    print STDERR "Expr $1\n";
    $filter = Maltools::Vcf::Expression->new(vcf => $vcf, expr => $1);
  \}
  return [$filter,$sign eq "+"];
\}
