package MalariaGEN::AGV::Job;

use strict;
use warnings;
use File::Temp;

use Carp qw(confess);
use Moose;

has 'wd' => ( is => 'ro', isa => 'Str', required => 1, 
    default => sub { $_ = `pwd`; chomp; $_; }, writer => '_set_wd');
has 'tool' => ( is => 'ro', isa => 'MalariaGEN::AGV::Tool', required => 1 );
has 'engine' => ( is => 'ro', isa => 'MalariaGEN::AGV::Engine', writer => '_set_engine' );
has 'inputs' => ( is => 'ro', isa => 'HashRef', default => sub { {} });
has 'outputs' => ( is => 'ro', isa => 'HashRef', default => sub {{}});
has 'cpu_ratio' => ( is => 'rw', isa => 'Num|Undef', default => 1);
has 'cpu_count' => ( is => 'rw', isa => 'Num|Undef', lazy => 1, default => 1);
has 'return_code' => ( is => 'rw', isa => 'Num|Undef', default => undef);
has 'error_message' => ( is => 'rw', isa => 'Str', default => "");
has 'tmp_required' => (is => 'rw', isa => 'Num|Undef', default => undef);
has 'name' => (is => 'rw', isa => 'Str', lazy => 1, builder => '_build_name');

# Expressed in Seconds:
# Wall time is then calculated as cpu_time / cpu_ratio 
has 'cpu_time' => ( is => 'rw', isa => 'Num|Undef', default => 60 * 60);
# Expressed in Megs:
has 'memory' => ( is => 'rw', isa => 'Num|Undef', default => 100); 
has 'is_distributed' => ( is => 'ro', isa => 'Bool', default => 0);

sub wall_time {
  return $_[0]->cpu_time / $_[0]->cpu_ratio;
}

sub has_input {
  my ($self,$name) = @_;
  return exists $self->inputs()->{$name};
}

sub data {
  my ($self,$name) = @_;
  my $tool_data = $self->tool->data(name => $name, as => 'REF');
  unless (defined $tool_data) {
    die "there is no tool data with name '$name'";
  }
  elsif ($tool_data->mode =~ /out/) {
    return $self->output($name);
  }
  else {
    return $self->input($name);
  }
}

sub input {
  my ($self,$name) = @_;
  unless ($self->has_input($name) ) {
    my $input = $self->tool()->input($name);
    confess "requested unknown input '$name'" unless $input; 
    return $self->_process_default($input->default());
  }

  my $result = $self->{inputs}->{$name};
  unless (defined($result)) {
     return wantarray ? () : undef;
  }
  elsif (ref($result) eq "ARRAY") {
     return wantarray ? @$result : $result; 
  }
  else {
     return wantarray ? ($result) : $result;
  }
}

sub set_data {
  my ($self,$data,$value) = @_;
  if (ref($data) eq "") {
     my $name = $data;
     $data = $self->tool->data(name => $name, as => 'REF') 
        or die "not such a tool data '$name'";
  }
  if ($data->mode eq 'in') {
    return $self->set_input($data->name => $value);
  } 
  else {
    return $self->set_output($data->name => $value);
  }
}

sub set_input {
  my ($self,$name,$value) = @_;
  $self->inputs()->{$name} = $value;
}





sub set_output {
  my ($self,$name,$value) = @_;
  $self->outputs()->{$name} = $value;
}

sub add_input {
  my ($self,$name,$value) = @_;
  $self->add_input_with_type($name,$value,'string');
}

sub uses_standard_io {
  my ($job) = @_;
  my $tool = $job->tool;
  my %file_data = map { ($_->name => $_ ) } grep { $_->type->has_files } values %{$tool->data()};
  foreach my $d (values %file_data) {
    my $v = $job->data($d->name);
    if (ref($v) eq "ARRAY") {
       foreach my $vv (@$v) { return 1 if defined $vv && $vv eq "-" };
    }
    else {
       return 1 if defined $v && $v eq "-";
    }
  }
  return 0;
}

sub add_input_with_type {
  my ($self,$name,$value,$type) = @_;
  my $type_name  = blessed($type) ?  $type->name : $type;
  if ($self->has_input($name)) {
    die "you cannot added an input with a name of an existing input";
  }
  $self->tool->add_input(name => $name, type => $type);
  $self->set_input($name,$value); 
}

sub output {
  my ($self,$name) = @_;
  return $self->outputs()->{$name} if ($self->has_output($name));
  die "requested non-existent output $name" unless defined $self->tool()->output($name);
  return $self->_process_default($self->tool()->output($name)->default());
}

sub has_output {
  my ($self,$name) = @_;
  return exists $self->outputs()->{$name};
}

sub command {
  my ($self,%args) = @_;
  $args{J} = \$self unless exists $args{J};
  $args{E} = \$self->engine unless exists $args{E};
  return $self->tool()->command(%args);
}

sub command_template {
  my ($self,%args) = @_;
  $args{J} = \$self unless exists $args{J};
  $args{E} = \$self->engine unless exists $args{E};
  return $self->tool()->command_template(%args);
}

sub script {
  my ($self,%args) = @_;
  $args{J} = \$self unless exists $args{J};
  $args{E} = \$self->engine unless exists $args{E};
  my $result = $self->tool()->script(%args);
  return $result;
}

sub script_template {
  my ($self,%args) = @_;
  $args{J} = \$self unless exists $args{J};
  $args{E} = \$self->engine unless exists $args{E};
  my $result = $self->tool()->script_template(%args);
  $result;
}

sub interpreter {
  my ($self,%args) = @_;
  $args{J} = \$self unless exists $args{J};
  $args{E} = \$self->engine unless exists $args{E};
  return $self->tool()->interpreter(%args);
}

sub _process_default {
  my ($self,$default) = @_;
  return wantarray ? () : undef unless defined $default;
  my $tool = $self->tool;
  my $engine = $self->engine;
  if (blessed($default) && $default->isa("Text::Template")) {
    return wantarray ? ($default->fill_in(HASH => { J => \$self, T => \$tool, E => \$engine})) : 
        $default->fill_in(HASH => { J => \$self, T => \$tool, E => \$engine}); 
  }
  elsif (ref($default) eq "") {
   return wantarray ? ($default) : $default;
  }
  elsif (ref($default) eq "CODE") {
   return $default->(J => $self, T => $tool, E => $engine);
  }
  else {
   return $default;
  }
}

sub tempfile {
 my ($self,$template,$suffix) = @_;
 $template ||= "tmp_XXXXX";
 $suffix ||= "";
 my ($tfh,$tfn) = File::Temp::tempfile($template, DIR => $self->wd);
 $tfh->close();
 unlink($tfn);
 return $tfn . $suffix;
}

sub tempdir {
 my ($self,$template) = @_;
 $template ||= "tmpdir_XXXXX";
 my ($tfn) = File::Temp::tempdir($template, DIR => $self->wd);
 return $tfn;
}

sub _build_name {
  my $self = shift;
  my $tool = $self->tool;
  my $class =ref($tool);
  $class =~ s/MalariaGEN::AGV::Tools:://;
  $class =~ s/::/\-/g;
  $class =~ s/(\b)([A-Z])/$1\u$2/g;
  return $class; 
}

1;
