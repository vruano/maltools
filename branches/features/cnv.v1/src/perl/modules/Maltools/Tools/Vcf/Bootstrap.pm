package Maltools::Tools::Vcf::Bootstrap;

use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
  in => { type => 'vcf', mandatory => 1},
  repeats => { type => 'num', mandatory => 1},
  size => { type => 'num', mandatory => 1},
  seed => { type => 'num', mandatory => 1},
  chr => { type => 'string', mandatory => 0},
  pattern => { type => 'string', mandatory => 0, default => 'bootstrap_%05d.vcf' },
};

our $OUTPUTS = {
  out_dir => {type => 'directory', mandatory => 1},
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

use Vcf;
use Fcntl qw(:seek);
use POSIX;
use IO::File;
use File::Spec::Functions qw(catfile);
use FileCache;
use Heap::Simple;

my $input = "{$J->input("in")}" || "-";
my $seed = {$J->input("seed") || "undef"};
my $repeats = {$J->input("repeats")};
my $size = {$J->input("size")};
my $out_dir = "{$J->output("out_dir")}";
my $pattern = "{$J->input("pattern")}";
my $summary_file = "{$J->output("summary")}" || undef;
my $chr = "{$J->input("chr")}" || undef;

my $fh = $input eq "-" ? \*STDIN : IO::File->new($input,'r') 
	or die "could not open input vcf file";

my $vcf = Vcf->new(fh => $fh);

$vcf->parse_header();

print STDERR "Calculating variant offsets (. = +16384) ";
my @offset = (tell($fh));
my $count = 0;
my $in_chr = defined $chr ? 0 : 1; 
while (my $line = $vcf->next_line) \{
 if (defined $chr) \{
   my $line_chr = substr($line,0,index($line,"\t"));
   if ($in_chr && $line_chr ne $chr) \{ # reached end of chromosome
     $in_chr = 0;
     last;
   \}
   $in_chr = $line_chr eq $chr;
   next unless $in_chr;
 \}
 $count++;
 push @offset,tell($fh); 
 print STDERR '.' unless $count & 16383; # 2^14  - 1;
\}
print STDERR "done\n";


srand($seed) if defined $seed;
my $count = scalar(@offset);

# flat size between 0 and 1 is interpreted 
# as the fraction of the total number of variants
if (int($size) != $size) \{
 if ($size > 1) \{
    $size = int($size);
 \}
 else \{
    $size = int($size * $count);
 \}
\}

print STDERR "Calculating repeat positions (. = +64)";
#Calculating positions per repeat.
my @repeats = (); $#repeats = $repeats - 1;
my @files = (); $#files = $repeats - 1;
my @already_opened = (); $#already_opened = $repeats - 1; 
for (my $i = 0; $i < $repeats; $i++) \{
  my $idx = floor(rand($count)); 
  my $offset = $offset[$i];
  my @positions = (); $#positions = $size;
  for (my $j = 0; $j < $size; $j++) \{
     $positions[$j] = int(rand($count)); 
  \}
  @positions = sort \{ $b <=> $a \} @positions;
  $repeats[$i] = \@positions;
  $files[$i] = catfile($out_dir,sprintf($pattern,$i));
  $already_opened[$i] = 0;
  print STDERR '.' unless $i & 63;
\}
print STDERR "done\n";

my $heap = Heap::Simple->new(order => '<', elements => 'Any');
for (my $i = 0; $i < $repeats; $i++) \{
  my $positions = $repeats[$i];
  $heap->key_insert($$positions[$#$positions],$i); 
\}

print STDERR "Generating bootstrap files (. ~= 10%) ";
my $remaining = $size * $repeats;
my $total = $remaining;
my $ten_pc = int($total / 10);
my $ten_approx = 1;
while ($ten_approx < $ten_pc) \{
  $ten_approx *= 2;
\}
$ten_pc = $ten_approx;
#Reading in variants and printing out bootstrap files.
while ($heap->count > 0) \{
  my $top_repeat = $heap->extract_top;
  my $positions = $repeats[$top_repeat]; 
  my $next_pos = pop @$positions;
  $heap->key_insert($$positions[$#$positions],$top_repeat) if ($#$positions >= 0);
  my $offset = $offset[$next_pos];
  my $ao = $already_opened[$top_repeat]++;
  my $mode = $ao ? '>>' : '>';
  my $out_fh = cacheout ($mode,$files[$top_repeat]);
  print $out_fh $vcf->format_header unless $ao;
  seek($fh,$offset,SEEK_SET);
  my $line = <$fh>;
  print $out_fh $line;
  $remaining--;
  print STDERR "."  unless $remaining & ($ten_pc - 1);
\}
print STDERR "done\n";

#Closing bootstrap files
foreach my $f (@files) \{
  cacheout_close($f);
\}

if ($summary_file) \{
  print STDERR "Generating summary file '$summary_file'\n";
  require JSON::XS;
  my $sfh = IO::File->new($summary_file,"w") or die "could not create summary file '$summary_file'";
  my $summary = \{
    repeats => $repeats,
    size => $size,
    seed => $seed,
    variants_count => $count,
    files => [@files]
  \};
  my $coder = JSON::XS->new->ascii->pretty->allow_nonref;
  print $sfh $coder->encode($summary);
  $sfh->close();
\}
else \{
  print STDERR "No summary file generated\n";
\}

