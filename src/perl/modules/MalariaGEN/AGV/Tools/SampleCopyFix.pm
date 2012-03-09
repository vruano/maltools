package MalariaGEN::AGV::Tools::SampleCopyFix;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

our $INPUTS = {
   "in" => { type => { name => 'bam', indexed => 1 },  mandatory => 1},
};

our $OUTPUTS = { 
    out => { type => { name => 'bam', indexed => 1 }, mandatory => 1 },
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 500;
}

sub calculate_cpu_time {
  return 10 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return "perl";
}

1;
__DATA__
{  $in   = $J->input('in');
  $out  = $J->output('out'); '' }
use strict;
use warnings;


my ($in_sample,$out_sample) = ("{$in}","{$out}");

my $new_header_file = $out_sample . ".hdr";
my $old_header = `samtools view -H $in_sample`;
my @old_header_lines = split(/\n/,$old_header);

my $new_header = join("\t",qw(@HD VN:1.3 SO:coordinate)) . "\n"; 

$new_header .= join("\n", grep \{ index($_,'@PG') < 0 \} @old_header_lines) . "\n";
$new_header .= join("\t",qw(@PG ID:bwa)) . "\n";

open(HDR,">$new_header_file");
print HDR $new_header;
close(HDR);

`samtools reheader $new_header_file $in_sample > $out_sample`;
if ($?) \{
  die "fail to fix $in_sample into $out_sample";
\}
unlink($new_header_file);

`samtools index $out_sample`;
if ($?) \{
  die "fail to index $out_sample";
\}

