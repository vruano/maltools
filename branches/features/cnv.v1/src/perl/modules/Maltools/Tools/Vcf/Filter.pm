package Maltools::Tools::Vcf::Filter;

use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  op => { type => 'string', mandatory => 1, default => 'add' },
  filters => {type => 'string', mandatory => 1},
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
use Maltools::Vcf::Filter;

my $input = "{ $J->input("in") }";
my $output = "{ $J->output("out") }";
my @filter_specs = qw({ join(" ",@\{$J->input("filters")\})});

open(STDIN,$input) or die "cannot open input file $input" if $input ne "-";
open(STDOUT,">$output") or die "cannot open output file $output" if $output ne "-";


my $vcf = Vcf->new(fh => \*STDIN);
$vcf->recalc_ac_an(0);
$vcf->parse_header();

my @filters = map \{ Maltools::Vcf::Filter->instance($_) \} @filter_specs;

foreach my $filter (@filters) \{
  $filter->add_header_line($vcf);
\}

print $vcf->format_header();

while (my $x = $vcf->next_data_hash) \{
  my $filters = join(";",@\{$x->\{FILTER\}\});
  my %filters_editions = ();
  foreach my $filter (@filters) \{
    $filters_editions\{$filter->id\} = 1 unless $filter->passes($x); 
  \}
  my $new_filters = [split(/;\s*/,$vcf->add_filter($filters,%filters_editions))];
  $x->\{FILTER\} = $new_filters;
  print $vcf->format_line($x);
\}
$vcf->close();

close(STDIN) if $input ne "-";
close(STDOUT) if $output ne "-";

