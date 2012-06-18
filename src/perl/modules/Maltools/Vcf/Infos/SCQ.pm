package Maltools::Vcf::Infos::SCQ;

use strict;
use warnings;

use Maltools::Vcf::CoverageAnalyzer;
use Moose;
use File::Spec::Functions qw(catfile);

extends 'Maltools::Vcf::Info';

has '+id' => (default => 'SCQ');
has '+description' => (lazy => 1, default => sub { 'Indicates sample coverage at particular genotype quality thresholds, concretelly: ' . $_[0]->thr });
#has '+number' => ( default => '.');
has '+type' => ( default => 'Float');
has 'thr' => (is => 'ro', isa => 'Str|ArrayRef', default => sub { '10,30,100,300' });
has 'gs_file' => (is => 'ro', isa => 'Str|Undef', default => undef );
has '_threshold_array' => (is => 'ro', isa => 'ArrayRef', lazy => 1 , default => sub { my $thr = $_[0]->thr; $thr = [split(/\,/,$thr)] if ref($thr) eq ""; return $thr; });
has '_counter_groups' => (is => 'rw', isa => 'HashRef', init_arg => undef, default => sub {{}});
has '_total_sample_count' => (is => 'rw', isa => 'Num|Undef', init_arg => undef, default => undef);
has 'good_samples_hash' => (is => 'rw', isa => 'HashRef', lazy => 1, builder => '_good_samples_hash_builder' );
has '_vcf' => (is => 'rw', isa => 'Maybe[Vcf]');

sub init_travesal {
  my ($self,$vcf) = @_;
  $self->SUPER::init_travesal($vcf);
  $self->_vcf($vcf);
  $self->_total_sample_count(scalar(keys %{$self->good_samples_hash}));
}

sub calculate_values {
  my ($self,$h,$vcf) = @_;
  my $gtypes = $h->{gtypes};
  my @thresholds = @{$self->_threshold_array};
  my @passed = map { 0 } @thresholds;
  my $ts = $self->_total_sample_count;
  my $gs = $self->good_samples_hash;
  foreach my $k (keys %$gtypes) {
    next unless exists($gs->{$k});
    my $v = $gtypes->{$k};
    my $gq = $v->{GQ} or next;
    for (my $i = 0; $i <= $#thresholds; $i++) {
      $passed[$i]++ if $gq > $thresholds[$i];
    }
  }
  #@passed = map { sprintf("%1.4f",$_/$ts) } @passed; 
  return [@passed];
}

sub _good_samples_hash_builder {
  my $self = shift;
  my $gs_file = $self->gs_file;
  my $vcf = $self->_vcf or die "good samples builder cannot be called before starting the travesal";
  my @samples = $vcf->get_samples;
  my %samples = map { $_ => 1 } @samples;
  my %included = ();
  if ($gs_file) {
     my $gs_fh = IO::File->new($gs_file,'r') or die "could not find the good sample file\n";
     while (my $line = <$gs_fh>) {
       next if $line =~ /^#/;
       next unless $line =~ /^\s*(\S+)*/;
       my $id = $1;
       exists($samples{$id}) or die "cannot find sample with name '$id' in input\n";
       $included{$id} = 1;
     }
     $gs_fh->close();
  }

  if (scalar(keys %included) == 0) {
     return { map { $_ => 1 } @samples };
  }
  else {
     return \%included;
  }

}



1;
