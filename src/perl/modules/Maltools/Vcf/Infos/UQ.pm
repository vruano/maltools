package Maltools::Vcf::Infos::UQ;

use strict;
use warnings;

use Maltools::Vcf::CoverageAnalyzer;
use Moose;
use File::Spec::Functions qw(catfile);

extends 'Maltools::Vcf::Info';

has '+id' => (default => 'UQ');
has '+description' => (default => 'Uniqueness score');
has '+number' => ( default => 1);
has '+type' => (default => 'Integer');
has 'uq_file' => (is => 'ro', isa => 'Str', required => 1);

has '_uq_file_fh' => (is => 'rw', isa => 'FileHandle', init_arg => undef, lazy => 1 , builder => '_uq_file_fh_builder');
has '_last_uq_line' => (is => 'rw', isa => 'Str', init_arg => undef, lazy => 1, builder => '_last_uq_line_builder');

sub _uq_file_fh_builder {
  my ($self) = @_;
  my $file_name = $self->uq_file;
  -e $file_name or die "cannot find uniqueness file $file_name";
  return IO::File->new($file_name,"r");
}

sub _last_uq_line_builder {
  my ($self) = @_;
  my $fh = $self->_uq_file_fh;
  my $line = <$fh>;
  return $line;
}

sub _read_uq_line {
  my $fh = $_[0]->_uq_file_fh;
  my $line = <$fh>;
  $_[0]->_last_uq_line($line);
  return $line;
}

sub calculate_values {
  my ($self,$h) = @_;
  my $seq = $h->{CHROM};
  my $pos = $h->{POS};
  return $self->uniqueness_score($seq,$pos);
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

sub uniqueness_score {
  my ($self,$seq,$pos) = @_;
  my $line = $self->_last_uq_line;
  while ($line) {
    chomp $line;
    my ($l_seq,$l_pos,$score) = split(/\t/,$line);
    my $cmp = $self->compare_seq_pos($seq,$pos,$l_seq,$l_pos);
    if ($cmp == 0) {
      return $score;
    }
    elsif ($cmp < 0) {
      return undef;
    }
    else {
      $line = $self->_read_uq_line;
    }
  }
  return undef;
}


1;
