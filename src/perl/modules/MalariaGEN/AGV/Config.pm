=pod

=head1 NAME

MalariaGEN::AGV::Config -- wired configuration file for the AGV pipeline

=head1 DESCRIPTION

This file contains all the basic stateless configuration for the AGV pipeline. 
It exports the configuration using convenience methods and configuration subclasses.

E<13>

=head2 Internals

The configuration data resides in a HashRef called $config. Each entry represents the configuration
of a different component or interface to external component. Each of these can be accessed using the 
object returned by the associated "config" method. E.g. for configuration details about our interface
with the Sanger's solexa pipeline:

use MalariaGen::AGV::Config qw(sanger_config);

my $sc = sanger_config();

my $run_file_location = $sc->lane_file(run => 5031,  lane_file => '5031_2_2_nonhuman');
# $run_file_location = '/fuse/mpsafs/runs/5031/5031_2_2_nonhuman.fastq'

The following list enumerates and  describes the configuration sections one by one

=over 4

=item sanger
E<13>
WTSI solexa pipeline interface configuration.

=item reference
E<13>
Reference genome sequence configuration

=item jobs
E<13>
Jobs execution settings including execution engine and program configurations.

=back

=head1 TODO

Externalize the configuration data to a text file, 
xml document or a database for an easier customization 
and language-independent quering.

=cut
package MalariaGEN::AGV::Config;


use warnings;
use strict;
use base 'Exporter';
use String::Interpolate qw(interpolate);
use JSON::XS;
use IO::File;
use File::Spec::Functions qw(catfile);
use Cwd qw(realpath getcwd);
use File::Basename qw(basename dirname);
use MalariaGEN::AGV::Manifest::Config;

our @EXPORT_OK = qw(scratch_config sanger_config reference_config jobs_config data_config config_reader config);

our $VERSION = '0.01';

my $AGV_HOME = $ENV{AGV_HOME} || $ENV{PGV_HOME} || guess_home(); 

sub template_dir {
  my ($class,@others) = @_;
  return resolve_config()->get_resource('templates',@others);
}

our $config = resolve_config();

sub new {
   return $config;
}


sub resolve_config {
  my $program = basename($0);
  
  my $home = $AGV_HOME || guess_home();

  my @files = (
    catfile(getcwd(),"CONFIG.json"),
    catfile(getcwd(),$program . ".json"),
    catfile($ENV{HOME},".pgv",$program . ".json"),
    catfile($home,"conf",$program. ".json"),
    catfile($ENV{HOME},".pgv","default.json"),
    catfile($home,"conf","default.json")
  );

  my $conf_file = undef;

  my $softdir = $ENV{PGV_HOME} || ".";
  my $result = MalariaGEN::AGV::Manifest::Config->new();
  $result->set_variable(program => $program, softdir => $softdir);
  
  foreach my $f (@files) {
    next unless -f $f;
    $result->load(file => $f);
  }
  
  return $result;
}




sub data_config {
  return DataConf->new({root => $config->get_datadir()});
}

sub scratch_config {
  return ScratchConf->new(data => $config->get('execution/scratch'));
}

sub reference_config {
  return ReferenceConf->new({config => $config});
}

sub jobs_config {
  return JobsConf->new($config->{jobs});
}


{
package DataConf;
use File::Spec::Functions qw(catfile);

sub new {
  my ($class,$data) = @_;
  my %data = %$data;
  my $self = \%data;
  bless $self,$class;
}

sub data_file {
  my ($self,@parts) = @_;
  my $file = catfile(@parts);
  return $file if index($file,$self->{root}) == 0;
  return catfile($self->{root},$file);
}

sub data_root {
  my ($self,@parts) = @_;
  return $self->{root};
}

sub map_loc {
  my ($self,%args) = @_;
  return catfile("maps",$args{run},join("_",$args{run},$args{lane}) . ".bam");
}

sub sample_dir {
  return raw_sample_dir(@_);
}

sub raw_sample_dir {
  my $self = shift;
  return catfile($self->{root},"samples","bam");
}

sub lane_dir {
  my $self = shift;
  return catfile($self->{root},'lanes');
}

sub lane_bam {
  my $self = shift;
  my $lane_name;
  if ($#_ == 0) {
    $lane_name = shift;
  }
  else {
    $lane_name = $_{run} .  "."  . $_{lane};
  }
  $lane_name =~ s/\./_/g;
  return catfile($self->lane_dir,$lane_name . ".bam")
}

sub aligned_sample_dir { 
  my $self = shift;
  return catfile($self->{root},"samples","aligned");
} 

sub realigned_sample_dir {
  my $self = shift;
  return catfile($self->{root},"samples","realigned");
}

sub recalibrated_sample_dir {
  my $self = shift;
  return catfile($self->{root},"samples","recalibrated");
}

sub sample_bam {
  my ($self,%args) = @_;
  my $result = catfile($self->{root},"samples",$args{class} || 'realigned',$args{ox_code} . ".bam");
  -f $result or return undef;
  if ($args{indexed}) {
    my $bai = $result . ".bai";
    -f $bai or return undef;
  }
  return $result;
}

sub genotyping_out_dir {
  my ($self,%args) = @_;
  my $result = catfile($self->{root},"genotypes",$args{group},$args{class} || 'realigned');
  return $result;
}

1;
}

{
package ScratchConf;
use Moose;

has 'data' => (is => 'ro', isa => 'HashRef', default => sub { {} });

sub scratch_instance {
  my $self = shift;
  my $data = $self->data;
  my $areas = $data->{areas};
  my @areas = ();
  for (my $i = 0; $i <= $#$areas; $i += 2) {
     my $type = $$areas[$i];
     my $root = $$areas[$i+1];
     $type = uc(substr($type,0,1)) . substr($type,1);
     my $class = "MalariaGEN::AGV::Engines::Scratch::$type";
     eval "require $class" or
         die "could not resolve scratch area type '$type'";
     push @areas, $class->new(root => $root);
  }
  require MalariaGEN::AGV::Engines::Scratch::Array
      or die "could not resolve scratch 'array' type";
  return MalariaGEN::AGV::Engines::Scratch::Array->new(areas => \@areas, exists($data->{exclusions}) ? ( exclusions => $data->{exclusions}) : ());
}

1;
}

{
package ReferenceConf;
use File::Spec::Functions qw(catfile);

sub new {
  my ($class,$data) = @_;
  my %data = %$data;
  my $self = \%data;
  bless $self,$class;
}

sub file_name {
  my ($self,%args) = @_;
  return $self->{config}->get_reference_sequences();
}

sub coding_tables {
  my ($self,%args) = @_;
  my $version = $args{version};
  return $self->{config}->_get_reference_satellite('.coding-tables');
}

sub uniqueness_tables {
  my ($self,%args) = @_;
  return $self->{config}->_get_reference_satellite('.uniqueness-tables');
}

sub possible_snps_list {
  my $self = shift;
  return $self->{config}->get_reference(@_)->possible_snps_file();
}

sub known_snps_mask_file {
  my $self = shift;
  return $self->{config}->get_reference(@_)->probable_snps_file();
}

sub known_snps_mask_format {
  my ($self) = @_;
  return $self->{config}->get_reference('knownSnpsMask','format');
}

1;
}

{
package JobsConf;
use Scalar::Util qw(blessed);
use MalariaGEN::AGV::Engine;
sub new {
  my ($class,$data) = @_;
  my %data = %$data;
  my $self = \%data;
  bless $self,$class;
}

sub engine {
  my $self = shift;
  my $engine = $self->{config}->get_default_execution_engine();
  return $engine unless ref($engine) eq "";
  return MalariaGEN::AGV::Engine->engine_by_name($engine);
}


1;
}

sub guess_home {
  return realpath(catfile(dirname($0),"..")); 
}


1;
