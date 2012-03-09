package MalariaGEN::AGV::Pipelines::Genvar;

use Moose;

use File::Basename qw(dirname);
use File::Path qw(make_path);

extends 'MalariaGEN::AGV::Pipeline';

override 'create' => sub {
  my $self = shift;
  my %args = @_;
   
  my @return = super();

  my $manifest = $self->manifest;
  $self->check_manifest($manifest);
  make_path $manifest->get_path('alignments');
  make_path $manifest->get_path('coverage');
  make_path $manifest->get_path('snpsdata');

  return @return;
};

sub check_manifest {
  my $self = shift;
  my $manifest = shift;

  my $aln_path = $manifest->get_path('alignments');
  $aln_path or die "there is no alignment path specified (empty 'paths/alignments') in the manifest!!!";
  $aln_path = dirname($aln_path); 
  -d $aln_path or die "the alignment parent path specified '$aln_path' is not a directory !!!";

  my $cvg_path = $manifest->get_path('coverage');
  $cvg_path or die "there is no coverage path specified (empty 'paths/coverage') in the manifest!!!";
  $cvg_path = dirname($cvg_path);
  -d $cvg_path or die "the coverage parent path specified '$cvg_path' is not a directory !!!";
  
  my $snps_path = $manifest->get_path('snpsdata');
  $snps_path or die "there is no snpsdata path specified (empty 'paths/snpsdata') in the manifest!!!";
  $snps_path = dirname($snps_path);
  -d $snps_path or die "the snpsdata parent path specified '$snps_path' is not a directory !!!";
  
  my $candidate_snps  = $manifest->get_parameter('candidateSnps');
  my $candidate_snp_list = $manifest->get_parameter('canidateSnpList');

  if ($candidate_snps ne 'all' && ! $candidate_snp_list ) {
     die "you must provide a candidate snp list file name unless 'parameters/candidateSnps' == 'all'";
  }

  # samples
  return 1;
  my $samples = $manifest->get_all_samples();

  if (ref($samples) eq "") {
    die "samples specification string or url '$samples' was not resolved in the manifest";
  }
  elsif (ref($samples) eq "ARRAY") {
    my @samples = grep { $_ } @$samples;
    if (scalar(@$samples) == 0) {
      die "there is no samples with data specified in the manifest";
    }
    else {
      my @scalar_samples = grep { ref($_) eq "" } @samples;
      die "there are samples was data has not been resolved in the manifest: ". join (", ",@scalar_samples);
    }
  }
}

1;
