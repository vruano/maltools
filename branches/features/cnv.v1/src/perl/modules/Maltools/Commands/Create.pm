package Maltools::Commands::Create;

use Moose;

extends 'Maltools::Command';

use Cwd qw(realpath);
use File::Basename qw(basename);
use Getopt::Long qw(GetOptionsFromArray);
use Maltools::Manifest::Pipeline;
use Maltools::Pipeline;

sub help_summary {
   return 'creates different data or pipelines';
}

sub help_text {
   my $self = shift;
   my $cl_name = $self->cl_name;
   my $summary = $self->help_summary;
   return <<EOM
Command:

  create - $summary

Synopsis:

  $cl_name create <object-type> options

  where <object-type> can be one of the following:

    'pipeline'.

Pipelines:

  $cl_name create pipeline -m manifest [-c class] -b basedir

      creates a pipeline of a given type (might be indicated within the manifest).

  -m json-file  indicates the configuration of the pipeline (inputs, outputs, 
                task to run, etc...). The content is dependent of the 
                pipeline type that might be indicated within he manifest 
                or explicilty using the -t option.

  -c string     (optional if provided within the manifest), indicates the pipeline
                class/type thus how to inteprete its content (e.g. how to complete
                its contents, what templates to use to create the execution bed).

  -b directory  the basedir, where to place the pipeline files
                of the pipeline instance.
EOM
;
}

sub execute {
  my ($self,$site,$params) = @_;
  my @args = @{$params->{arguments}};
  unless (@args) {
     return $self->error_return("you must specify and object type as a first argument");
  }
  my $object_type = lc(shift @args);
  if ($object_type eq 'pipeline') {
     return $self->pipeline(@args);
  }
  else {
     return $self->error_return("unknown object type '$object_type'");
  }
}

sub pipeline {
  my $self = shift;
  my $manifest = undef;
  my $type = undef;
  my $dest_dir = undef;
  my %vars = ();
  GetOptionsFromArray(\@_,"manifest|m=s" => \$manifest,"class|c=s" => \$type,"basedir|b=s" => \$dest_dir, "variable|v=s%" => \%vars);

  return $self->error_return('you must provide a manifest file using option -m') unless ($manifest);
  return $self->error_return("manifest file '$manifest' does not exists or is not accesible") unless -e $manifest;
  return $self->error_return("manifest file '$manifest' is not a regular file") unless -f $manifest;
  return $self->error_return("manifest file '$manifest' cannot be read, wrong permissions") unless -r $manifest;

  return $self->error_return('you must provide a destination directory using option -b') unless $dest_dir;
  return $self->error_return("destination directory '$dest_dir' is in fact not a directory") if -e $dest_dir && ! -d $dest_dir;
 
  my $pipeline_class; 
  if ($type) {
    $pipeline_class = Maltools::Pipeline->class_by_name($type)
        or return $self->error_return("unknown pipeline type '$type'"); 
  }

  my $manifest_object = Maltools::Manifest::Pipeline->new(file => $manifest,
                       variables => { %vars, softdir => $ENV{MALTOOLS_HOME}, outdir => '.', basedir => $dest_dir, prog => basename($0)} )
    or return $self->error_return("could not load manifest '$manifest'");

  unless($type) {
    $type = $manifest_object->get('class')
        or return $self->error_return("the manifest '$manifest' does not indicate the pipeline class (/class element), you need to indicate one explicitly");
    $pipeline_class = Maltools::Pipeline->class_by_name($type)
        or return $self->error_return("unknown pipeline type '$type' indicated in the manifest");
  }

  my $pipeline = $pipeline_class->new(manifest => $manifest_object); 
  $pipeline->create(basedir => $dest_dir, %vars);
  return $self->ok_return();
}

1;
