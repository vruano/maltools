package MalariaGEN::AGV::Pipelines::Probgen;

use Moose;
extends 'MalariaGEN::AGV::Pipeline';

use File::Basename qw(dirname);
use File::Path qw(make_path);

has '+manifest' => (is => 'ro', isa => 'MalariaGEN::AGV::Pipelines::Probgen::Manifest', required => 1);

override 'create' => sub {
  my $self = shift;
  my %args = @_;
   
  my @return = super();

  my $manifest = $self->manifest;
  $self->check_manifest($manifest);
  make_path $manifest->get_path('alignments');


  return @return;
};

sub check_manifest {
  my $self = shift;
  my $manifest = shift;

  my $aln_path = $manifest->get_path('alignments');
  $aln_path or die "there is no alignment path specified (empty 'paths/alignments') in the manifest!!!";
  $aln_path = dirname($aln_path); 
  -d $aln_path or die "the alignment parent path specified '$aln_path' is not a directory !!!";

  return 1;
}


1;