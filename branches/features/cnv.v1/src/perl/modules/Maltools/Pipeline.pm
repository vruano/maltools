package Maltools::Pipeline;

use Moose;

use File::Find qw(find);

use Maltools::Config;
use Maltools::Manifest::Pipeline;

use IO::File;
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use File::Basename qw(basename dirname);
use Template;
use Cwd qw(getcwd);

has manifest => (is => 'ro', isa => 'Maltools::Manifest::Pipeline', required => 1);

around BUILDARGS => sub {
  my $orig = shift;
  my $class = shift;
  my %args = @_;
  my $manifest = $args{manifest};
  return $class->$orig(@_) unless $manifest || blessed($manifest);
  $args{manifest} = _BUILDARGS_process_manifest($class,$manifest);
  return $class->$orig(%args);
};

sub _BUILDARGS_process_manifest {
  my $class = shift;
  my $manifest = shift;
  return undef unless defined $manifest;
  if (ref($manifest) eq "") {
    return $class->_manifest_class->new(file => $manifest);
  }
  elsif (ref($manifest) eq "HASH") {
    return $class->_manifest_class->new(%$manifest);
  }
  elsif (ref($manifest) eq "ARRAY") {
    return $class->_manifest_class->new(@$manifest);
  }
  elsif (blessed($manifest) && $manifest->isa('Maltools::Manifest')) {
    if ($manifest->isa($class->_manifest_class)) {
      return $manifest;
    }
    else {
      return $class->_manifest_class->meta->rebless_instance($manifest->clone());
    }
  }
  else {
    die "cannot handle manifest reftype '" . ref($manifest) . "'";
  }
}

sub class_by_name {
  my $class = shift;
  my $name = shift or return undef;

  if ($name =~ /^Maltools::Pipeline/) {
    eval "require $name" and return $name;
    return undef;
  }
  elsif ($name =~ /^::/) {
    $name = 'Maltools::Pipelines' . $name;
    eval "require $name" and return $name;
    return undef;
  }
  elsif ($name =~ /::/) {
    $name = 'Maltools::Pipelines::' . $name;
    eval "require $name" and return $name;
    return undef;
  }
  else {
    $name = 'Maltools::Pipelines::' . uc(substr($name,0,1)) . substr($name,1);
    eval "require $name" and return $name;
    return undef;
  }
}


sub class {
  my $self = shift;
  my @parts = split(/::/,ref($self));
  my $last_name = $parts[$#parts];
  return lc($last_name);
}

sub _manifest_class {
  my $class = shift;
  print STDERR $class,"\n";
  my $manifest_attr = $class->meta->get_attribute('manifest');
  unless ($manifest_attr->has_type_constraint) {
    return 'Maltools::Manifest::Pileline';
  }
  else {
    my $candidate = $manifest_attr->type_constraint->name;
    return 'Maltools::Manifest::Pipeline' unless $candidate;
    
    eval "require $candidate"
      or die "cannot load the manifest type '$candidate'";
    return $candidate;
  }
}

sub _template_dir {
   my $result = Maltools::Config->template_dir("pipeline",$_[0]->class);
   return $result; 
}

sub _template {
 my $self = shift;
 my $template = Template->new({
     DEBUG => 1,
     EVAL_PERL => 1,
     ABSOLUTE => 1,
     INCLUDE_PATH => $self->_template_dir(),
 });
 return $template;
}

sub create {
  my $self = shift;
  my %args = @_;
  my $dest_dir = $args{basedir};

  make_path $dest_dir;

  my $manifest = $self->manifest->clone();

  $manifest->auto_complete(%args);

  $manifest->save(catfile($dest_dir,"MANIFEST.json"));

  $self->apply_templates(manifest => $manifest, destdir => $dest_dir);
}

sub run {
  my $self = shift;
  my %args = @_;
  my $basedir = $args{basedir}
     or die "no base-directory provided";
  -d $basedir or die "the bese-directory does not seem to exists"; 
  my $cur_dir = getcwd();
  chdir $basedir;
  eval { $self->_run(%args) };
  my $error = $@;
  chdir $cur_dir;
  if ($error) {
    die $error;
  }
}

sub _run {
  my $self = shift;
  my %args = @_;
  my $basedir = $args{basedir} || getcwd();
  my $out_file = catfile($basedir,'output.txt');
  my $err_file = catfile($basedir,'error.txt');
  my $threads = $args{cpus} || 200;
  my $make_opts = $args{mkopts} || '-k';
  my $target = $args{target} || 'all';
  print STDERR "make -C $basedir -j $threads $make_opts $target 2> $err_file 1> $out_file\n";
  system("make -C $basedir -j $threads $make_opts $target 2> $err_file 1> $out_file");
  $? and die "Errors executing the pipeline (code=$?)";
}

sub apply_templates {
  my $self = shift;
  my %args = @_;
  my $dest_dir = $args{destdir};
  my $manifest = $args{manifest} || $self->manifest;

  my $template = $self->_template();

  my $template_dir = $self->_template_dir();

  my $vars = {
     manifest => $manifest,
     basedir => $dest_dir,
     program => basename($0),
  };
  find({
    no_chdir => 1,
    wanted => sub {
       my $fn = $File::Find::name;
       my $rel_fn = $fn;
       $rel_fn =~ s/\Q$template_dir\E\/*//;
       return if ($fn =~ /^\./ || $fn =~ /\.inc$/);
       if (-d $fn) {
         make_path(catfile($dest_dir,$rel_fn));
       }
       elsif (-f $fn) {
         my $dest_fn = catfile($dest_dir,$rel_fn); 
         $dest_fn =~ s/\.tt2?$//;
         make_path(dirname($dest_fn));
         open(my $dest_fh,">",$dest_fn);
         $template->process($fn,$vars,$dest_fh)
             or die "Error in template: " . $template->error;
         $dest_fh->close();
         my $perm = (stat $fn)[2] & 07777;
         chmod($perm, $dest_fh);
         close($dest_fn);
       }
    },
  },$template_dir);

}  



1;
