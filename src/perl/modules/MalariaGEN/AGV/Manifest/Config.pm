package MalariaGEN::AGV::Manifest::Config;

use Moose;
use MalariaGEN::AGV::Reference;

extends 'MalariaGEN::AGV::Manifest';

sub get_reference {
  my $self = shift;
  my $ref_spec = $self->get('reference',@_);
  if (ref($ref_spec) eq "HASH") {
    $ref_spec = $self->_json2perl_key_translate($ref_spec);
  }
  elsif (ref($ref_spec) eq "ARRAY") {
    $ref_spec = $self->_json2perl_key_translate({@$ref_spec});
  }
  return MalariaGEN::AGV::Reference->new($ref_spec);
}

sub get_reference_sequences {
  my $self = shift;
  my $ref = $self->get_reference(@_);
  return $ref->sequence_file;
}

sub get_datadir {
  my $self = shift;
  return $self->get('data') || $self->get('variables/datadir');
}

sub get_default_execution_engine {
  my $self = shift;
  return $self->get('execution/defaultEngine') || 'local';
}


sub get_execution_engine {
  my $self = shift;
  my $engine = shift;
  return $self->get('execution/engines',$engine);
}


sub get_reference_annotation {
  my $self = shift;
  return $self->get_reference(@_)->annotation_file;
}

1;
