package MalariaGEN::AGV::Config::Reader;

use Moose;

use JSON::XS;
use IO::File;
use Storable;
use String::Interpolate qw(interpolate);

has 'config' => (is => 'ro', isa => 'HashRef', default => sub { {} });

sub load {
   my $self = shift;
   my $new_config = $self->_read_file(shift);
   $self->_merge_config($new_config);
}

sub clone {
   my $self = shift;
   $self->meta->name->new(config => Storable::clone($self->config())); 
}

sub vars {
  return $_[0]->get('vars') || {};
}

sub set_var {
  my $self = shift;
  my $key = shift;
  ref($key) && die "the variable ref type must be an scalar";
  my $value = shift || "";
  ref($value) && die "the variable value type must be an scalar";
  my $vars = $self->get('vars');
  if (!$vars) {
     $vars = $self->get()->{vars} = {};
  }
  elsif (ref($vars) ne "HASH") {
     die "illegal vars ref-type " . ref($vars);
  }
  $vars->{$key} = $value;
}

sub set_vars {
  my $self = shift;
  my %new_vars = @_;
  foreach my $k (keys %new_vars) {
     my $v = $new_vars{$k};
     $self->set_var($k,$v);
  } 
}

sub get {
  my $self = shift;
  my $path = join("/",@_);
  my @path_elements = split(/\//,$path);
  my $return = $self->config();

  foreach my $pe (@path_elements) {
    if (ref($return) ne "HASH") {
       return undef;
    }
    $return = $return->{$pe};
    return undef unless defined $return;
  }
  return $return;
}

sub _read_file {
  my $self = shift;
  my $file_name = shift;
  my $json = JSON::XS->new();

  my $fh = IO::File->new($file_name,'r') or die "could not read '$file_name'";

  while (my $line = <$fh>) {
     $json->incr_parse($line);
  }

  my $obj = $json->incr_parse;

  my $variables = $obj->{vars} || $obj->{variables} || {};

  ref($variables) eq "HASH" or die "the variables reference type must be a hash";
  my $interpolables = { %{$self->vars()} };
  my $new_variables = {};
  foreach my $k (keys %$variables) {
     my $v = $$variables{$k};
     my $new_v = interpolate($v,$interpolables);
     $$new_variables{$k} = $new_v;
     $$interpolables{$k} = $new_v;
  }
  $variables = $new_variables;

  foreach my $k (keys %$obj) {
     next if $k eq 'vars';
     $$obj{$k} = _process_config_element($interpolables,$$obj{$k});
  }

  $$obj{vars} = $variables;

  return $obj;
}

sub _merge_config {
  my $self = shift;
  my $incoming = shift;
  my $current = shift || $self->config();
  foreach my $k (keys %$incoming) {
     my $in_value = $$incoming{$k};
     my $cu_value = $$current{$k};
     unless (defined $cu_value) {
       $$current{$k} = $in_value;
       next;
     }    
     my $in_rt = ref($in_value);
     my $cu_rt = ref($cu_value);
     if ($in_rt ne $cu_rt || !$in_rt) {
       $$current{$k} = $in_value;
       next;
     }
     elsif ($in_rt eq "ARRAY") {
       push @{$$current{$k}},@$in_value;
     }
     elsif ($in_rt eq "HASH") {
       $self->_merge_config($in_value,$cu_value);
     }
     else {
       die "don't know how to handle ref-type '$in_rt'";
     }
  }
}

sub _process_config_element {
  my $vars = shift;
  my $element = shift;
  my $ref_type = ref($element);

  if ($ref_type eq "ARRAY") {
    return _process_config_array($vars,$element);
  }
  elsif ($ref_type eq "HASH") {
    return _process_config_hash($vars,$element);
  }
  elsif (!$ref_type) {
     return defined $element ? interpolate($element,$vars) : undef;
  }
  else {
     die "cannot handle ref-type '$ref_type'";
  }

}

sub _process_config_hash {
  my $vars = shift;
  my $hash = shift;
  my $new_hash = {};
  foreach my $k (keys %$hash) {
     my $v = $$hash{$k};
     my $new_k = interpolate($k,$vars);
     $$new_hash{$new_k} = _process_config_element($vars,$v);
  }
  return $new_hash;
}

sub _process_config_array {
  my $vars = shift;
  my $array = shift;
  $array = [ map { _process_config_element($vars,$_) } @$array ];
  return $array;
}


1;
