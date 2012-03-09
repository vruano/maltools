package MalariaGEN::AGV::Vcf::Infos::CODING;

use strict;
use warnings;

use Moose;
use File::Spec::Functions qw(catfile);
use MalariaGEN::AGV::GffUtils::Iterator;

extends 'MalariaGEN::AGV::Vcf::Info';

has '+id' => (default => 'CODING');
has '+description' => (default => 'Flags coding sites');
has '+number' => ( default => undef);
has '+type' => (default => 'Flag');
has 'gff_file' => (is => 'ro', isa => 'Str', required => 1);

has '_gff_iterator' => (is => 'ro', init_arg => undef, lazy => 1 , builder => '_gff_iterator_builder');

sub _gff_iterator_builder {
  my ($self) = @_;
  my $file_name = $self->gff_file;
  -e $file_name or die "cannot find Gff file $file_name " . `pwd; ls -l $file_name`;
  return MalariaGEN::AGV::GffUtils::Iterator->new(file => $file_name);
}

sub calculate_values {
  my ($self,$h) = @_;
  my $seq = $h->{CHROM};
  my $pos = $h->{POS};
  my $gff_iterator = $self->_gff_iterator;
  $gff_iterator->seek_pos($seq,$pos); 
  my $features = $gff_iterator->current_features;
  foreach my $f (@$features) {
    my $type = $f->{type};
    return 1 if $type eq "exon";
    return 1 if $type eq "CDS";
    return 1 if $type eq "cds";
  }
  return 0;
}

sub compare_seq_pos {
  my ($self,$seq,$pos,$oseq,$opos) = @_;
  
  if ($seq eq $oseq) {
    return $pos <=> $opos; 
  }
  elsif ($seq eq "MT") {
    return 1;
  }
  elsif ($oseq eq "MT") {
    return -1;
  }
  else {
    return $seq cmp $oseq;
  }
}

1;
