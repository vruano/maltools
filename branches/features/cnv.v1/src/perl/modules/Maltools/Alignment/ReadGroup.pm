package Maltools::Alignment::ReadGroup;

use Moose;

our $NEXT_ID = "RG_0000";

has id => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_id' );
has line => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_line', predicate => '_has_line');
has sequencing_center => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_sequencing_center');
has description => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_description');
has date => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_date');
has flow_order => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_flow_order');
has key_sequence => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_key_sequence');
has library => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_library');
has platform => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_platform');
has platform_unit => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_platform_unit');
has sample => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_sample');
has alignment => (is => 'ro', isa => 'Maltools::Alignment', required => 1);


sub forward_fastq_file {
  return $_[0]->alignment->_forward_fastq_file($_[0]);
}

sub reverse_fastq_file {
  return $_[0]->alignment->_reverse_fastq_file($_[0]);
}

sub _build_line {
  my $self = shift;
  my @parts = ();
  
  push @parts, _line_part('ID',$self->id);
  push @parts, _line_part('CN',$self->sequencing_center);
  push @parts, _line_part('DS',$self->description);
  push @parts, _line_part('DT',$self->date);
  push @parts, _line_part('FO',$self->flow_order);
  push @parts, _line_part('KS',$self->key_sequence);
  push @parts, _line_part('LB',$self->library);
  push @parts, _line_part('PL',$self->platform);
  push @parts, _line_part('PU',$self->platform_unit);
  push @parts, _line_part('SM',$self->sample);
} 

sub _line_part {
   my $tag = shift;
   my $value = shift;
   return $value ? ("${tag}:${value}") : ();
}

sub _build_from_line {
  my $self = shift;
  my $tag = shift;
  
  if ($self->_has_line) {
    my $line = $self->line;
    return $1 if $line =~ /\t${tag}:([^\t]+)/;
  }
  return undef;

}

sub _build_id {
  return $_[0]->_build_from_line('ID') || $NEXT_ID++;
}

sub _build_sequencing_center {
  return $_[0]->_build_from_line('CN');
}

sub _build_description {
  return $_[0]->_build_from_line('DS');
}

sub _build_date {
  return $_[0]->_build_from_line('DT');
}

sub _build_flow_order {
  return $_[0]->_build_from_line('FO');
}

sub _build_key_sequence {
  return $_[0]->_build_from_line('KS');
}

sub _build_library {
  return $_[0]->_build_from_line('LB');
}

sub _build_platform {
  return $_[0]->_build_from_line('PL');
}

sub _build_platform_unit {
  return $_[0]->_build_from_line('PU');
}

sub _build_sample {
  return $_[0]->_build_from_line('SM');
}


1;
