package Maltools::Vcf::Info;

use strict;
use warnings;
use Moose;

use Vcf;

has 'id' => (is => 'ro', isa => 'Str', required => 1);
has 'number' => (is => 'ro', isa => 'Maybe[Int]', default => undef);
has 'type' => (is => 'ro', isa => 'Str', default => 'String');
has 'description' => (is => 'ro' , isa => 'Str', default => undef);

sub init_travesal {
  my ($self,$vcf) = @_;
}

sub finish_travesal {
  my ($self,$vcf) = @_;
}

sub add_header_line {
  my ($self,$vcf) = @_;
  return $vcf->add_header_line({key => 'INFO', ID => $self->id, Number => defined($self->number) ? $self->number : -1, 
      Type => $self->type, Description => $self->description});
}

sub calculate_values {
  my ($self,$x) = @_;
  die "Not implemented yet";
}

sub process_site {
  my ($self,$x) = @_;
  my $values = $self->calculate_values($x);
  return 0 unless defined $values;
  if ($self->type eq 'Flag') {
    $x->{INFO}->{$self->id} = undef if $values;
  }
  elsif(ref($values) eq "ARRAY") {
    $x->{INFO}->{$self->id} = join(",",@$values);
  }
  else {
    $x->{INFO}->{$self->id} = $values;
  }
  return 1;
}

sub info_class_by_id {
  my ($class,$id) = @_;
  my $cls = "Maltools::Vcf::Infos::$id";
  eval "require $cls" or die "cannot find info with id $id";
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
  my $cls = $class->info_class_by_id($id);
  return $cls->new(%spec);
}

sub _instance_from_scalar {
  my ($class,$spec) = @_;
  die "cannot resolve filter from undefined spec" unless defined ($spec); 
  my ($id,@params)  = split(/:/,$spec);
  
  my %params = map { ( split(/\s*=\s*/,$_) ) } @params;
  my $cls = $class->info_class_by_id($id);
  return $cls->new(id => $id, %params);
}

1;
