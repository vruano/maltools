package MalariaGEN::AGV::Engines::Local;

use strict;
use warnings;

use Moose;
extends 'MalariaGEN::AGV::Engine';
use MalariaGEN::AGV::Tool;
use Cwd qw(getcwd realpath);
use File::Temp;
use MalariaGEN::AGV::Config qw(sanger_config);
use File::Spec::Functions qw(catfile file_name_is_absolute);
use IO::File;
use File::Basename qw(dirname);
use File::Path qw(remove_tree make_path);

has '+name' => ( default => 'Local');
has 'root' => ( is => 'ro', isa => 'Str' , lazy => 1, default => sub { getcwd });
has 'temp_root' => (is => 'ro', isa => 'Str', lazy => 1, default => sub { my $result = catfile($_[0]->root,"tmp"); make_path($result); return $result; });
has '_sc' => ( is => 'ro' , lazy => 1, default => sub { sanger_config });

sub new_script_file {
  my ($self,$cwd) = @_;
  return $self->tempfile(template=>'script_XXXX');
}

sub run_job {
  my ($self,$job,%args) = @_; 
  my $tool = $job->tool;
  $job->_set_engine($self);
  my $procs = $args{procs} || 1;
  my $cwd = $self->tempdir(dir => $self->root);
  $job->_set_wd($cwd);
  my $sfn = $self->new_script_file($cwd);
  my $sfh;
  $self->_prepare_inputs($job);
  $self->_prepare_outputs($job);
  my %job_args = (
     J => \$job, T => \$tool, E => \$self, D => $cwd,
     S => \$sfn, MAKE => "lsmake" , SH => "tcsh", P => $procs, %args); 
  my $script = $job->script(%job_args);
  $sfh = IO::File->new($sfn,"w");
  print $sfh $script;
  $sfh->close();
  my $cmd = $job->command(%job_args);
  if ($cmd =~ /^\$MAKE/) {
    $cmd =~ s/\$MAKE/make/;
    $cmd .= " -C $cwd";
  }
  else {
    $cmd =~ s/^\$SHELL/tcsh/;
  }
  eval { print STDERR $cmd,"\n...";system($cmd) };
  print STDERR "\n";
  my $code = $?;
  my $error = $@;
  eval { remove_tree($cwd); } unless $code || $@;
  $job->return_code($code);
  $job->error_message($@);
  return $code == 0 && !$error;
}

sub tempdir {
  my ($self,%args) = @_;
  my $template = $args{template} || "plain_XXXXXXX";
  delete $args{template};
  return File::Temp::tempdir($template, DIR => $self->temp_root,%args);
}

sub tempfile {
  my ($self,%args) = @_;
  my $unlink = $args{unlink} || 0;
  delete $args{unlink};
  my $template = $args{template} || "plain_XXXXXXX";
  delete $args{template};
  my ($fh,$fn) = File::Temp::tempfile($template, DIR => $self->temp_root,%args);
  $fh->close();
  unlink $fn if $unlink;
  return $fn;
}

sub _prepare_inputs {
  my ($self,$job) = @_;
  my $tool = $job->tool;
  my %file_inputs = map { ($_->name => $_ ) } grep { $_->type->has_files } values %{$tool->inputs()};
  my $uses_standard_io = 0;
  foreach my $in (values(%file_inputs)) {
    my $value = $job->input($in->name) or next;
    my @files = $in->type->file_names($value);
    my $svalue = undef;
    foreach my $ofile (@files) {
      $uses_standard_io = 1 if $ofile eq "-";
      next if $ofile eq "-";
      my $file = realpath($ofile) || $ofile;
      $svalue = $svalue || $in->type->relocate_file_on_value($value,$ofile,$file);
    }
    $job->set_input($in->name => $svalue) if defined $svalue;
  }
  return $uses_standard_io;
}

sub _prepare_outputs {
  my ($self,$job) = @_;
  my $tool = $job->tool;
  my %file_inputs = map { ($_->name => $_ ) } grep { $_->type->has_files } values %{$tool->outputs()};
  my $uses_standard_io = 0;
  foreach my $in (values(%file_inputs)) {
    my $value = $job->output($in->name) or next;
    my @files = $in->type->file_names($value);
    my $svalue = undef;
    foreach my $ofile (@files) {
      $uses_standard_io = 1 if $ofile eq "-";
      next if $ofile eq "-";
      my $file = realpath($ofile) || $ofile;
      my $dirname = dirname($file);
      -e $dirname or make_path($dirname);
      $svalue = $svalue || $in->type->relocate_file_on_value($value,$ofile,$file);
    }
    $job->set_output($in->name => $svalue) if defined $svalue;
  }
  return $uses_standard_io;
}


1;
