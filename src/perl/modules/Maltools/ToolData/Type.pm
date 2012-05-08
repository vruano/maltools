package Maltools::ToolData::Type;

use File::Spec::Functions qw(catfile splitpath catpath);
use Moose;

has 'name' => ( is => 'ro', isa => 'Str');
has 'has_files' => ( is => 'ro', default => 0 );

sub type_class_by_name {
  my ($class,$name) = @_;
  my $cls = "Maltools::ToolData::$name";
  eval "require $cls";
  $@ and die "cannot find type tool-data type $name: $@";
  return $cls;
}

sub extends_type {
  my ($self,$candidate) = @_;
  return $self->isa("Maltools::ToolData::$candidate");
}

sub instance {
  my ($class,%spec) = @_;
  my $type = $spec{name} or die "type spec does not contain name";
  my $cls = $class->type_class_by_name($spec{name});
  my $cls_meta = $cls->meta();
  my %cls_specs = %spec;
  foreach my $k (keys %cls_specs) {
    delete $cls_specs{$k} unless $cls_meta->get_attribute($k);
  }
  return $cls->new(%spec);
}

sub _relocate_array_value {
  my ($self,$value,$prev,$next) = @_;
  return $value if $prev eq $next;
  return [ map { return $self->relocate_value ($_,$prev,$next) } @$value ];
}

sub relocate_file_on_value {
  my ($self,$value,$prev,$next) = @_;
  if (ref($value) eq "") {
    return $self->relocate_value($value,$prev,$next);
  }
  else { # assume is ARRAY 
    for (my $i = 0; $i <= $#$value; $i++) {
       eval { $$value[$i] = $self->relocate_value($$value[$i],$prev,$next) };
    }
    return $value;
  }
}

sub relocate_value {
  my ($self,$value,$prev,$next) = @_;
  if (ref($value) eq "ARRAY") {
     return $self->_relocate_array_value($value,$prev,$next);
  } 
  return $value if $prev eq $next;

  my @prev = splitpath($prev);
  my @value = splitpath($value);
  my @next = splitpath($next);
  
  while ($#prev >= 0 && $#value >=0 && $prev[0] eq $value[0]) {
    shift @prev; shift @value;
  }
  while ($#next >= 0 && $#prev >=0 && $prev[$#prev] eq $next[$#next]) {
    pop @next; pop @prev;
  }
  return catpath(@next,@value) if $#prev == -1;
  die "cannot relocate '$value' based on '$prev' to '$next'";
}

1;
