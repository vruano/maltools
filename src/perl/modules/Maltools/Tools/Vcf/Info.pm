package Maltools::Tools::Vcf::Info;

use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  op => { type => 'string', mandatory => 1, default => 'add' },
  skip_filtered => { type => 'bool', default => 1 },
  infos => {type => 'string', mandatory => 1},
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
use Maltools::Vcf::Info;

my $input = "{ $J->input("in") }";
my $output = "{ $J->output("out") }";
my $skip_filtered = { $J->input("skip_filtered") };
my @info_specs = qw({ join(" ",@\{$J->input("infos")\})});

open(STDIN,$input) or die "cannot open input file $input" if $input ne "-";
open(STDOUT,">$output") or die "cannot open output file $output" if $output ne "-";


my $vcf = Vcf->new(fh => \*STDIN);
$vcf->recalc_ac_an(0);
$vcf->parse_header();

my @infos = map \{ Maltools::Vcf::Info->instance($_) \} @info_specs;

foreach my $info (@infos) \{
  $info->add_header_line($vcf);
\}
print $vcf->format_header();


foreach my $info (@infos) \{
  $info->init_travesal($vcf);
\}
while (my $x = $vcf->next_data_hash) \{
  my $process_it = 1;
  if ($skip_filtered) \{
    my $filters = uc(join(";",@\{$x->\{FILTER\}\}));
    $process_it = 0 unless $filters eq "." || $filters eq "" || $filters eq "PASS";
  \}
  if ($process_it) \{ 
    foreach my $info (@infos) \{
      eval \{ $info->process_site($x); \};
      $@ and die $info->id," ",$@,"\n";
    \}
  \}
  print $vcf->format_line($x);
\}

foreach my $info (@infos) \{
  $info->finish_travesal($vcf);
\}

$vcf->close();

close(STDIN) if $input ne "-";
close(STDOUT) if $output ne "-";

