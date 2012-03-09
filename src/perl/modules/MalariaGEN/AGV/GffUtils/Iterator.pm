package MalariaGEN::AGV::GffUtils::Iterator;

use Moose;
use strict;
use warnings;
use IO::File;

has 'fh' => (is => 'ro', isa => 'FileHandle', lazy => 1, default => sub {{ defined($_[0]->file) ? IO::File->new($_[0]->file) : undef }});
has 'file' => (is => 'ro', isa => 'Str'); 
has 'current_pos' => (is => 'ro', isa => 'Num', writer => '_set_current_pos');
has 'current_seq' => (is => 'ro', isa => 'Str', writer => '_set_current_seq');
has '_current_features' => (is => 'rw', isa => 'ArrayRef' , 
             lazy => 1, init_arg => undef, default => sub { [] } );
	has '_last_line' => (is => 'rw', isa => 'Str|Undef', lazy => 1,
             init_arg => undef, builder => '_read_next_line' );

sub seek_pos {
  my ($self,$seq,$pos) = @_;
  $seq or die "a sequence name must be provided as fist argument";
  $pos or die "a position must be provided as second argument";
  if (defined $self->current_pos) {
    return 1 if $self->current_pos eq $pos && $self->current_seq eq $seq;
  }
  $self->_read_to_pos($seq,$pos);
}

sub current_features {
  my ($self,$type) = @_;
  my $candidate = $self->_current_features;
  my @result = @$candidate;
  if (defined $type) {
    @result = grep { $_->{type} eq $type } @result;
  }
  return wantarray ? @result : \@result;
}

sub _read_to_pos {
  my ($self,$seq,$pos) = @_;
  my $line = $self->_last_line;
  my @features = $self->current_features;
  @features = grep { 
      $_->{seq} eq $seq && $_->{start} <= $pos && $_->{stop} >= $pos } @features;
  my $in_seq = 0;
  while (defined $line) {
    if ($line !~ /^#/) {
      my ($l_seq,$db,$type,$start,$stop) = split(/\t/,$line);
      $in_seq ||= $seq eq $l_seq;
      last if ($l_seq ne $seq && $in_seq);
      last if ($start > $pos);
      if ($start <= $pos && $stop >= $pos) {
        push @features, { seq => $seq, start => $start, stop => $stop, type => $type, db => $db, line => $line  };
      }
    } 
    $line = $self->_read_next_line;
  }
  $self->_update_seq_and_pos($seq,$pos,\@features);
  $self->_current_features(\@features);
  return scalar(\@features);
}

sub _update_seq_and_pos {
  my ($self,$seq,$pos,$features) = @_;
  if (scalar(@$features) > 0) {
     $self->_set_current_seq($seq);
     $self->_set_current_pos($pos);
  }
  else {
     $self->_set_current_pos(-1);
     $self->_set_current_seq("");
  }
}

sub _read_next_line {
  my ($self) = @_;
  my $fh = $self->fh;
  my $nl = <$fh>;
  $self->_last_line($nl);
  return $nl;
}

1;


