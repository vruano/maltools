package Maltools::Vcf::Infos::CCP;

use strict;
use warnings;

use Maltools::Vcf::CoverageAnalyzer;
use Moose;
use File::Spec::Functions qw(catfile);

extends 'Maltools::Vcf::Info';

has '+id' => (default => 'CPP');
has '+description' => (default => 'Indicates the cumulative probability of the observed coverage in that position based on the site class (SC)');
has '+number' => ( default => 1);
has '+type' => ( default => 'Float');
has 'cvg_dir' => (is => 'ro', isa => 'Str', required => 1);
has '_counter_groups' => (is => 'rw', isa => 'HashRef', init_arg => undef, default => sub {{}});

our $ca = 'Maltools::Vcf::CoverageAnalyzer';

sub counter_group {
   my ($self,$seq) = @_;
   my $cgs = $self->_counter_groups;
   return $cgs->{$seq} if exists($cgs->{$seq});
   my $cvg_dir = $self->cvg_dir;
   my $file = catfile($cvg_dir,"cvg-$seq.json");
   -e $file or $file = catfile($cvg_dir,"cvg-all.json");
   my $fh = IO::File->new($file,"r");
   my $decoder = JSON::XS->new();
   while (my $line = <$fh>) {
    $decoder->incr_parse($line);
   }
   $fh->close();   
   return $cgs->{$seq} = $decoder->incr_parse();
}

sub calculate_values {
  my ($self,$h) = @_;
  my $info = $h->{INFO};
  my $gtypes = $h->{gtypes};
  my $total_ad = 0;
  foreach my $s (keys %$gtypes) {
    my $sample_ad = $gtypes->{$s}->{AD} || 0;
    my @counts = split(/,/,$sample_ad);
    foreach my $c (@counts) { $total_ad += $c };
  }
  my $coding = exists($info->{CODING});
  my $counter_group = $self->counter_group($h->{CHROM});
  my $counter = $counter_group->{$coding ? 'exonic' : 'non_coding'};
  
  my $result = $counter->{cumulative}->[$total_ad];
  print STDERR "$total_ad $result\n";
  return $result;
}


1;
