package MalariaGEN::AGV::ToolData;

use strict;
use warnings;
use Scalar::Util qw(blessed);
use MalariaGEN::AGV::ToolData::Type;
use Moose;
use Carp qw(longmess);

use Moose::Util::TypeConstraints qw(enum);

has 'name' => ( is => 'ro', isa => 'Str', required => 1 );
has 'type' => ( is => 'ro' , required => 1); 
has 'mandatory' => ( is => 'ro', isa => 'Bool', default => 0 );
has 'values' => ( is => 'ro', isa => 'ArrayRef' );
has 'default' => ( is => 'ro' );
has 'mode' => ( is => 'ro', isa => 'Str', required => 1 );


around BUILDARGS => sub {
   my ($orig,$class,%args) = @_;
   $args{type} = $class->_process_type($args{type},\%args);
   $class->$orig(%args);
};

sub _process_type {
  my ($class,$type,$args) = @_;

  return undef unless defined $type;
  if (ref($type) eq "") {
    return  MalariaGEN::AGV::ToolData::Type->instance(name => $type);
  }
  elsif (ref($type) eq "HASH") {
    return MalariaGEN::AGV::ToolData::Type->instance(%$type);
  }
  elsif (blessed($type)) {
    return $type if $type->isa('MalariaGEN::AGV::ToolData::Type');
  }
  die "cannot handle ref-type " . ref($type) . " as a tool-data type"; 
}

1;
