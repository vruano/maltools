package MalariaGEN::AGV::Tools::Vcf::Concat;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::DataTemplateTool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  inputs => { type => 'vcf', mandatory => 1, multiple => 1, mandatory => 1},
  reference => { type => 'reference', mandatory => 1 },
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

use IO::File;
use Getopt::Long;
use MalariaGEN::AGV::Reference;

my $reference = "{ $J->input('reference') }";
my @in_files = ({ "'". join("','", @{$J->input('inputs')}) . "'" });
my $out_file = "{ $J->output('out') }";

sub process_file_list \{
  my ($opt,$file_name) = @_;

  open my $fh , $file_name or die "could not open '$file_name' to read";
  my @inputs = <$fh>;
  close $fh;
  chomp @inputs;
  push @in_files, @inputs;
\}


my @bad_inputs = grep \{ not $_ && -f $_ && -r $_ \} @in_files;

die "invalid input files (not found, regular nor readable):\n\t" . join("\n\t",@bad_inputs) . "\n" if @bad_inputs;

die "must indicate a reference" unless defined $reference;
die "could not open to read or is not a regular file '$reference'" unless -f $reference && -r $reference; 

$reference = MalariaGEN::AGV::Reference->new(sequence_file => $reference);

my %positions = map \{ ($_ => peek_position ($_)) \} @in_files;

@in_files = sort \{ $positions\{$a\}->compare_with($positions\{$b\}) \} @in_files;
print STDERR join("\n",@in_files),"\n";
my $out_fh = defined $out_file && $out_file ne '-' ? IO::File->new($out_file,"w") : \*STDOUT;
$out_fh or die "could not open output $out_file\n";
print_vcf(shift(@in_files) => $out_fh, skip_header => 0);
print_vcf($_ => $out_fh, skip_header => 1) foreach @in_files;
close $out_fh;
exit 0;

sub print_vcf \{
  process_vcf($_[0],out_fh => $_[1],@_[2 .. $#_]);
\}

sub peek_position \{
  process_vcf($_[0],peek_start => 1);
\}

sub process_vcf \{
  my $file = shift;
  my %opts = @_;

  my $skip_header = $opts\{skip_header\} || $opts\{peek_start\};
  
  open my $in_fh , $file or die "could not read from $file\n";
  die unless defined $in_fh;
  my $line = <$in_fh>;
  if ($skip_header) \{
    while ($line) \{
      last unless index($line,'#') == 0;
  die unless defined $in_fh;
      $line = <$in_fh>;
    \}
  \}
  if ($opts\{peek_start\}) \{
    my ($chr,$pos) = split(/\s+/,$line);
    close $in_fh;
    return $line ? $reference->location($chr => $pos) : undef;
  \}
  my $out_fh = $opts\{out_fh\};
  while ($line) \{
    print $out_fh $line;
  die unless defined $in_fh;
    $line = <$in_fh>;
  \}
  close $in_fh;
\}




