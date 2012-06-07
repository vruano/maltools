package Maltools::Manifest;

use Moose;

use IO::File;
use JSON::XS;
use String::Interpolate;
use Storable;

has 'class' => (is => 'ro', isa => 'Str', lazy => 1, 
                predicate => 'has_class', writer => '_set_class', builder => '_build_class' );
has 'content' => (is => 'ro', isa => 'HashRef', lazy => 1, builder => '_build_content', predicate => '_has_content');
has 'file' => (is => 'ro', isa => 'Str', predicate => 'has_file');
has 'variable_hash' => (is => 'ro', isa => 'HashRef', lazy => 1, init_arg => 'variables', 
                   builder => '_build_variable_hash', writer => '_set_variable_hash');

sub BUILD {
    if ($_[0]->has_file) {
      $_[0]->load(file => $_[0]->file);
    }
    $_[0]->_class_update;
}

sub _class_update {
  my $self = shift;
  if (!$self->has_class) {
    my $class = $self->get("class");
    if ($class && !ref($class)) {
        $self->_set_class($class);
    }
  }
}

sub _build_class {
  my $self = shift;
  return $self->get("class") || 'unknown'; 
}

sub load {
   my $self = shift;
   my %options = @_;
   if ($options{file}) {
      $self->_merge_content($self->_read_file($options{file})) if exists $options{file};
   }
   if ($options{content}) {
      $self->_merge_content($options{content});
   }
   $self->_class_update;
}


sub save {
   my $self = shift;
   my $file = shift
       or die "you need to specify a file name";
   my $str = JSON::XS->new->pretty->encode($self->content);  
   my $fh = IO::File->new($file,'w')
       or die "could not open '$file' to write";
   print $fh $str;
   $fh->close;
}

our %_NARROWS = (
  pipe => 'Maltools::Manifest::Pipeline',
  pipeline => 'Maltools::Manifest::Pipeline'
);

sub narrow {
  my $self = shift;
  my $class = shift;
  unless($class && $self->has_class) {
    $class = $self->class;
  }
  $class = $_NARROWS{lc($class)} if exists $_NARROWS{lc($class)};    
  $class or die "you need to spicify a narrow class";
  eval "require $class" or die "there is no valid narrow manifest class '$class'";
  $class->isa('Maltools::Manifest')
    or die "'$class' is not decendant of Manifest base class";
  return $class->meta->rebless_instance($self->clone());
}


sub clone {
   my $self = shift;
   my %args = ();
   my $new_content = Storable::dclone($self->content());
   return $self->meta->clone_object($self,content => $new_content);
}


sub _build_variable_hash {
  my $self = shift;
  if ($self->_has_content) {
    my $result = $self->get('variables');
    unless (defined $result) {
        $result = {};
       $self->get()->{variables} = $result;
    } 
    return $result;
  }
  else {
    return {};
  }
}

sub get_variable {
  my $self = shift;
  my @names = @_;
  my $variables = $self->get('variables') || {};
  my @result = map { $variables->{$_} } @names;
  if ($#names == 0 && !wantarray) {
    return $result[0];
  }
  else {
    return wantarray ? @result : \@result;
  }
}

sub _build_content {
  return {};
}

sub _hash_translate {
   my $self = shift;
   my $orig = shift;
   my $trans = shift;

   my $result = {};
   foreach my $k (keys (%$orig)) {
     if (exists($trans->{$k})) {
       $result->{$trans->{$k}} = $orig->{$k};
     }
     else {
       $result->{$k} = $orig->{$k};
     }
   }
   return $result;
}

sub _json2perl_key_translate {
    my $self = shift;
    my $orig = shift;
    my $trans = shift;
    
    my $result = {};
    foreach my $k (keys (%$orig)) {
      if ($k =~ /[A-Z]/) {
        my @parts = ();
        while ($k =~ /([a-zA-Z1-9][A-Z1-9]*[a-z1-9]*)/g) { push @parts, $1 };
        @parts = map { uc($_) } @parts;
        my $key = join("_",@parts);
        $$result{$key} = $$orig{$k};
      }
      else {
        $$result{$k} = $$orig{$k};
      }
    }
    return $result;
}

sub _set_single_variable {
  my $self = shift;
  my $key = shift;
  ref($key) && die "the variable ref type must be an scalar";
  my $value = shift || "";
  ref($value) && die "the variable value type must be an scalar";
  my $vars = $self->get('variables');
  if (!$vars) {
     $vars = $self->get()->{variables} = {};
  }
  elsif (ref($vars) ne "HASH") {
     die "illegal vars ref-type " . ref($vars);
  }
  $vars->{$key} = $value;
}

sub set_variable {
  my $self = shift;
  my %new_vars = @_;
  foreach my $k (keys %new_vars) {
     my $v = $new_vars{$k};
     $self->_set_single_variable($k,$v);
  } 
}

sub get {
  my $self = shift;
  my $path = join("/",@_);
  my @path_elements = split(/\//,$path);
  my $return = $self->content();
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
  my $json = JSON::XS->new()->canonical;

  my $fh = IO::File->new($file_name,'r') or die "could not read '$file_name'";

  while (my $line = <$fh>) {
     $json->incr_parse($line);
  }

  my $obj = $json->incr_parse;

  my $variables = $obj->{variables} || $obj->{variables} || {};

  ref($variables) eq "HASH" or die "the variables reference type must be a hash";
  my $interpolables = { %{$self->variable_hash()} };
  my $new_variables = { %$interpolables };
  my @revisit = ();
  foreach my $k (keys %$variables) {
     my $v = $$variables{$k};
     next if exists $$interpolables{$k};
     my $new_v = $self->interpolate($v,$interpolables);
     $$new_variables{$k} = $new_v;
     push @revisit , $k;
  }
  
  $interpolables =  { %$new_variables };
  foreach my $k (@revisit) {
     my $v = $$variables{$k};
     my $new_v = $self->interpolate($v,$interpolables);
     $$new_variables{$k} = $new_v;
     $$interpolables{$k} = $new_v;
  }
  
  $variables = $new_variables;
  
  $self->_set_variable_hash($variables);

  foreach my $k (keys %$obj) {
     next if $k eq 'variables';
     $$obj{$k} = $self->_process_content_element($$obj{$k});
  }

  $$obj{variables} = $variables;

  return $obj;
}

sub _merge_content {
  my $self = shift;
  my $incoming = shift;
  unless (ref($incoming) eq "HASH") {
    die "content must be a hash but it had ref-type '" . ref($incoming) . "'";
  }
  my $current = shift || $self->content();
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
       next;
     }
     elsif ($in_rt eq "ARRAY") {
       push @{$$current{$k}},@$in_value;
     }
     elsif ($in_rt eq "HASH") {
       $self->_merge_content($in_value,$cu_value);
     }
     else {
       die "don't know how to handle ref-type '$in_rt'";
     }
  }
}

sub _process_content_element {
  my $self = shift;
  my $vars = $self->variable_hash;
  my $element = shift;
  my $ref_type = ref($element);

  if ($ref_type eq "ARRAY") {
    return $self->_process_content_array($element);
  }
  elsif ($ref_type eq "HASH") {
    return $self->_process_content_hash($element);
  }
  elsif (!$ref_type) {
     return defined $element ? $self->interpolate($element,$vars) : undef;
  }
  else {
     die "cannot handle ref-type '$ref_type'";
  }

}

sub _process_content_hash {
  my $self = shift;
  my $vars = $self->variable_hash;
  my $hash = shift;
  my $new_hash = {};
  foreach my $k (keys %$hash) {
     my $v = $$hash{$k};
     my $new_k = $self->interpolate($k,$vars);
     $$new_hash{$new_k} = $self->_process_content_element($v);
  }
  return $new_hash;
}

sub _process_content_array {
  my $self = shift;
  my $vars = $self->variable_hash;
  my $array = shift;
  $array = [ map { $self->_process_content_element($_) } @$array ];
  return $array;
}

sub interpolate {
  my $self = shift;
  my $str = shift;
  my $start = substr($str,0,2);
  my $a_pos = index($start,'!');
  return substr($str,1) if ($a_pos == 0);
  $str = substr($str,1) if ($a_pos == 1 && $start eq "\\!");
  return String::Interpolate::interpolate($str,@_);
}




1;
