package MalariaGEN::AGV::ROD;

use Moose;
use File::Spec::Functions qw(catfile);


has 'handle' => (is => 'ro', required => 1);
has '_last_line' => (is => 'rw', isa => 'Str', init_arg => undef, lazy => 1, builder => '_last_line_builder');

sub _last_line_builder {
  my ($self) = @_;
  my $fh = $self->handle;
  my $line = <$fh>;
  return $line;
}

sub _read_line {
  my $fh = $_[0]->handle;
  my $line = <$fh>;
  $_[0]->_last_line($line);
  return $line;
}

sub _compare_seq_pos {
  my ($self,$seq,$pos,$oseq,$opos) = @_;
  
  if ($seq eq $oseq) {
    return $pos <=> $opos; 
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
