package MalariaGEN::AGV::Vcf::Infos::CLASS;

use strict;
use warnings;

use Moose;
use File::Spec::Functions qw(catfile);

extends 'MalariaGEN::AGV::Vcf::Info';

has '+id' => (default => 'CLASS');
has '+description' => (lazy => 1, default => sub { 'Marks variants as credible or typable based on several threshold   ' . $_[0]->thr });
has '+number' => ( default => '1');
has '+type' => ( default => 'String');
has 'vc_qthr' => (is => 'ro' , isa => 'Maybe[Num]', default => 30);
has 'vt_qthr_hetero' => (is => 'ro', isa => 'Maybe[Num]', lazy => 1, default => sub { $_[0]->vc_qthr });
has 'vt_qthr' => (is => 'ro' , isa => 'Maybe[Num]', default => 30);
has 'vt_cthr' => (is => 'ro' , isa => 'Maybe[Num]', default => undef);
has 'svc_cvg_thr' => (is => 'ro', isa => 'Num' , default => 0.5);
has 'thr' => (is => 'ro', isa => 'Str|ArrayRef', default => sub { '10,30,100,300' });
has 'gs_file' => (is => 'ro', isa => 'Maybe[Str]', default => undef );
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
  return undef if $self->vc_qthr > $h->{QUAL};
  my $ts = $self->_total_sample_count;
  my $gs = $self->good_samples_hash;
  my $count = 0;
  my $vt_qthr = $self->vt_qthr;
  my $vt_cthr = $self->vt_cthr;
  my $vt_qthr_hetero = $self->vt_qthr_hetero;
  foreach my $k (keys %$gtypes) {
    next unless exists($gs->{$k});
    my $v = $gtypes->{$k};
    my $gt = $v->{GT};
    next if $gt =~ /\..\./;
    my ($a1,$a2) = split (/\\|\||\//,$gt);
    my $gq = $v->{GQ} or next;
    my $dp = $v->{DP} or next;
    next if $vt_qthr && $gq < $vt_qthr; 
    next if $vt_cthr && $dp < $vt_cthr;
    next if $vt_qthr_hetero && $a1 == $a2 && $gq < $vt_qthr_hetero;
    $count++;
  }
  my $cvg = $count / scalar(keys %$gs);
  return 'Typable' if ($cvg >= $self->svc_cvg_thr);
  return 'Credible';
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
