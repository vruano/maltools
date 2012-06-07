package Maltools::Commands::Run;

use Moose;

extends 'Maltools::Command';

use strict;
use warnings;

use Getopt::Long qw(GetOptionsFromArray);
use File::Spec::Functions qw(catfile);
use Maltools::Manifest::Pipeline;
use Maltools::Pipeline;

sub help_summary {
   return 'runs executable components such as pipelines';
}

sub help_text {
   my $self = shift;
   my $cl_name = $self->cl_name;
   my $summary = $self->help_summary;
   return <<EOM
Command:

  run - $summary

Synopsis:

  $cl_name run <object-type> options

  where <object-type> can be one of the following:

    'pipeline'.

Pipelines:

  $cl_name run pipeline -b <pipeline-installation-dir>

      runs a pipeline in its basedir or installation directory

  --basedir|-b where the target pipeline was created.

EOM
;
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV = @{$params->{arguments}};
  my $object_type = lc(shift @ARGV);
  if ($object_type eq 'pipeline') {
     return $self->pipeline(@ARGV);
  }
  else {
     return $self->error_return("unknown object type '$object_type'");
  }
}

sub pipeline {
  my $self = shift;
  my $type = undef;
  my $dest_dir = undef;
  my $threads = $self->running_options->{procs};
  my %vars = ();
  GetOptionsFromArray(\@_,"basedir|b=s" => \$dest_dir,"threads|t=i" => \$threads);

  return $self->error_return('you must provide a base directory using option -b') unless $dest_dir;
  return $self->error_return("the basedir provided '$dest_dir' does not exists") unless -e $dest_dir;
  return $self->error_return("the basedir provided '$dest_dir' is in fact not a directory") if -e $dest_dir && ! -d $dest_dir;
  
  my $manifest_file = catfile($dest_dir,'MANIFEST.json');
  
  return $self->error_return("cannot find a manifest file in the basedir '$dest_dir'") unless -e $manifest_file;
  return $self->error_return("the manifest file '$manifest_file' isn't in fact a regular file'") unless -f $manifest_file;
  return $self->error_return("the manifest file '$manifest_file' cannot be read, wrong permissions") unless -r $manifest_file;
 
  my $manifest_object = Maltools::Manifest::Pipeline->new(file => $manifest_file)
    or return $self->error_return("could not load manifest '$manifest_file'");
  
  $type = $manifest_object->get('class')
        or return $self->error_return("the manifest '$manifest_file' does not indicate the pipeline class (/class element), you need to indicate one explicitly");
  my $pipeline_class = Maltools::Pipeline->class_by_name($type)
        or return $self->error_return("unknown pipeline type '$type' indicated in the manifest");

  my $pipeline = $pipeline_class->new(manifest => $manifest_object); 
  $pipeline->run(basedir => $dest_dir, cpus => $threads);
  return $self->ok_return();
}

1;
