package MalariaGEN::AGV::Vcf::Filter;

use strict;
use warnings;
use Moose;

use Vcf;

has 'id' => (is => 'ro', isa => 'Str', required => 1);
has 'description' => (is => 'ro' , isa => 'Str', default => undef);

sub add_header_line {
  my ($self,$vcf) = @_;
  return $vcf->add_header_line({key => 'FILTER', ID => $self->id, Description => $self->description});
}

sub passes {
  return 1;
}


sub filter_class_by_id {
  my ($class,$id) = @_;
  my $cls = "MalariaGEN::AGV::Vcf::Filters::$id";
  eval "require $cls" or die "cannot find filter with id $id";
  return $cls;
}

sub instance {
  my ($class,@spec) = @_;
  if (scalar(@spec) == 1) {
    return $class->_instance_from_hash(%{$spec[0]}) if ref($spec[0]) eq "HASH";
    return $class->_instance_from_hash(@{$spec[0]}) if ref($spec[0]) eq "ARRAY";
    return $class->_instance_from_scalar($spec[0]) if ref($spec[0]) eq "";
  }
  else {
    return $class->_instance_from_hash(@spec);
  }
}

sub _instance_from_hash {
  my ($class,%spec) = @_;
  my $id = $spec{id};
  my $cls = $class->filter_class_by_id($id);
  return $cls->new(%spec);
}

sub _instance_from_scalar {
  my ($class,$spec) = @_;
  die "cannot resolve filter from undefined spec" unless defined ($spec); 
  my ($id,@params)  = split(/:/,$spec);
  
  my %params = map { ( split(/\s*=\s*/,$_) ) } @params;
  my $cls = $class->filter_class_by_id($id);
  return $cls->new(id => $id, %params);
}

1;
