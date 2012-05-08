package Maltools::Tools::Vcf::Coverage;

use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  gff => {type => 'gff', mandatory => 1},
  skip_filters => { type => 'flag', default => 1},
};

our $OUTPUTS = {
  out => {type => 'directory', mandatory => 1},
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

use strict ;
use warnings ;
use Vcf;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use Maltools::Vcf::CoverageAnalyzer;

my $gff_file = "{ $J->input("gff") }";
my $input_file = "{ $J->input("in") }";
my $out_dir = "{ $J->output("out") }";
my $skip_filtered = {$J->input("skip_filtered")};

#####
### Main body

my $coverage = Maltools::Vcf::CoverageAnalyzer->new(gff => $gff_file, filter_all => $skip_filtered); 

my $vcf = $input_file eq "-" ? Vcf->new(fh => \*STDIN) : Vcf->new(file => $input_file);
$vcf->parse_header;

my $counters = $coverage->analyze($vcf);

foreach my $s (keys %$counters) \{
 my $sample_counters = $counters->\{$s\};
 my $outdir = $s eq "all" ? $out_dir : catfile($out_dir,$s);
 -e $outdir or make_path($outdir);
 foreach my $k (keys %$sample_counters) \{
   my $out_file_prefix = catfile($outdir,"cvg-$k");
   my $json_file = $out_file_prefix . ".json";
   my $image_file = $out_file_prefix . ".png";
   $coverage->write_counter_group($sample_counters->\{$k\},$json_file);
   $coverage->coverage_plot($sample_counters->\{$k\},$image_file);
 \}
\}

