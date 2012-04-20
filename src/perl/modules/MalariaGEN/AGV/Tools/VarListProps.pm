package MalariaGEN::AGV::Tools::VarListProps;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use vars qw(%ENV);
use File::Basename qw(dirname);
use Text::Template;
use IO::File;
use Cwd qw(realpath);
use JSON::XS;
use File::Spec::Functions qw(catfile file_name_is_absolute);
use POSIX;

our $INPUTS = {
   "manifest" => { type => { name => 'file' },  mandatory => 1},
   "annotation" => { type => { name => 'gff' }, mandatory => 0},
   "literature_snps"  => { type => 'file'},
   "reference" =>  { type => { name => 'fasta', indexed => 1 } },
   "uniqueness_scores" => { type => 'file' },
   "snp_properties" => { type => 'file' },
   "coding_regions" => { type => 'file' },
   "coverage_cutoffs" => { type => 'file' },
   "candidate_snps" => { type => 'file' },
   "snp_lists" => { type => 'file', multiple => 1 },
   "coverage_files" => { type => 'file', multiple => 1},
};


our $OUTPUTS = { 
   "out" => { type => 'file' }
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 4000;
}

sub calculate_cpu_time {
  return 4 * 60 * 60; # for now just applicable to Pfalciparum 
}

sub interpreter {
   my ($self) = @_;
   return 'perl';
}

sub job { 
  my $self = shift;
  my $job = $self->SUPER::job(@_);
  my $manifest = $job->input('manifest');
  my $json = JSON::XS->new();
  my $mfh = IO::File->new($manifest,'r') or die "could not open '$manifest'";
  while (my $line = <$mfh>) { $json->incr_parse($line) };
  my $mani_object = $json->incr_parse();
  $mfh->close();
  my $base_path = realpath($mani_object->{paths}->{base_path} ||  dirname($manifest));
  my $annotation = $job->input('annotation') || $self->_apply_base_path($base_path,$mani_object->{paths}->{annotation});
  my $reference = $job->input('reference') || $self->_apply_base_path($base_path,$mani_object->{paths}->{reference});
  my $uniqueness_scores = $job->input('uniqueness_scores') || $self->_apply_base_path($base_path,$mani_object->{paths}->{uniqueness_scores});
  my $snp_properties = $job->input('snp_properties') || $self->_apply_base_path($base_path,$mani_object->{paths}->{snp_properties});
  my $literature_snps = $job->input('literature_snps') || $self->_apply_base_path($base_path,$mani_object->{paths}->{literature_snps});
  my $coverage_cutoffs = $job->input('coverage_cutoffs') || $self->_apply_base_path($base_path,$mani_object->{paths}->{coverage_cutoffs});
  my $coding_regions = $job->input('coding_regions') || $self->_apply_base_path($base_path,$mani_object->{paths}->{coding_regions});
  my $snpdata_path = $self->_apply_base_path($base_path,$mani_object->{paths}->{snpdata_path} || "." );
  my $coverage_path = $self->_apply_base_path($base_path,$mani_object->{paths}->{coverage_path} || ".");
  my $candidate_snps = $job->input('candidate_snps') || $self->_apply_base_path($base_path,$mani_object->{paths}->{candidate_snps});
  my $samples = $mani_object->{samples};
  $job->set_input('base_path',$base_path); 
  $job->set_input('annotation',$annotation);
  $job->set_input('reference',$reference);
  $job->set_input('uniqueness_scores',$uniqueness_scores);
  $job->set_input('coding_regions',$coding_regions);
  $job->set_input('snp_properties',$snp_properties);
  $job->set_input('literature_snps',$literature_snps);
  $job->set_input('coding_regions',$coding_regions);
  $job->set_input('coverage_cutoffs',$coverage_cutoffs);
  $job->set_input('candidate_snps',$candidate_snps);
  $job->set_output('out',$coverage_cutoffs);
  $job->set_input('snp_lists',[map { realpath (catfile($snpdata_path,$_ . ".snps")) } @$samples]); 
  $job->set_input('coverage_files',[map { realpath (catfile($coverage_path,$_ . ".coverage")) } @$samples]); 
  return $job;
}

sub _apply_base_path {
  my ($self,$base,$target) = @_;
  if (file_name_is_absolute($target)) {
    return realpath($target);
  }
  else {
    return realpath(catfile($base,$target));
  }
}


1;

__DATA__
{ 
  use File::Spec::Functions qw(catfile);
  use File::Basename qw(basename);
  $original_manifest  = $J->input('manifest');
  $base_path = "."; 
  $annotation = $J->input('annotation');
  $literature_snps  = $J->input('literature_snps');
  $reference = $J->input('reference');
  $uniqueness_scores = $J->input('uniqueness_scores');
  $snp_properties = $J->input('snp_properties');
  $coverage_cutoffs = $J->output('out');
  $coding_regions = $J->input('coding_regions');
  $candidate_snps = $J->input('candidate_snps');
  $snp_lists = $J->input('snp_lists');
  $coverage_files = $J->input('coverage_files');
  $manifest  = $J->tempfile(); 
  $coverage_path = catfile($J->wd(),'coverage'); 
  mkdir $coverage_path;
  foreach my $file (@$coverage_files) {
     my $bn = basename($file);
     if (-l $file) {
       $dest = catfile($coverage_path,$bn);
       `cp $file $dest`;
       $? and die "failed to copy $file into $dest"; 
     }
     else {
       symlink $file,catfile($coverage_path,$bn) or die "could not creae symbolic link from $file";
     }
  }

'' }

use strict;
use warnings;
use JSON::XS;
use IO::File;
use Cwd qw(realpath);
use File::Spec::Functions qw(catfile);

my $json = JSON::XS->new()->pretty();

my $original_manifest = "{$original_manifest}";
my $manifest = "{$manifest}";
my $coverage_path = "{$coverage_path}";


my $mfh = IO::File->new($original_manifest,'r');
my $mani_content = "";
while (my $line = <$mfh>) \{ $mani_content .= $line \};
$mfh->close();
my $mani_object = $json->decode($mani_content);

$mani_object->\{paths\}->\{base_path\} = "{$base_path}";
$mani_object->\{paths\}->\{annotation\} = "{$annotation}";
$mani_object->\{paths\}->\{literature_snps\} = "{$literature_snps}";
$mani_object->\{paths\}->\{reference\} = "{$reference}";
$mani_object->\{paths\}->\{uniqueness_scores\} = "{$uniqueness_scores}";
$mani_object->\{paths\}->\{snp_properties\} = "{$snp_properties}";
$mani_object->\{paths\}->\{coverage_cutoffs\} = "{$coverage_cutoffs}";
$mani_object->\{paths\}->\{coding_regions\} = "{$coding_regions}";
$mani_object->\{paths\}->\{candidate_snps\} = "{$candidate_snps}";
$mani_object->\{paths\}->\{out_path\} = "{$base_path}";

my $snp_list_example = "{$snp_lists->[0]}";
$snp_list_example =~ s/\/[^\/]+\.snps//;
$mani_object->\{paths\}->\{snpdata_path\} = $snp_list_example;
$mani_object->\{paths\}->\{coverage_path\} = $coverage_path;

my $coverage_cutoffs = $mani_object->\{paths\}->\{coverage_cutoffs\};

$mfh = IO::File->new($manifest,'w');
print $mfh $json->encode($mani_object);
$mfh->close();

$manifest = realpath($manifest);

my $pgv_home = $ENV\{MALTOOLS_HOME\};
my $vfwd = catfile($pgv_home,'resources','pgv_pipeline','variation_filtering2');


print STDERR "the coverage_cutoffs file is $coverage_cutoffs\n";
unlink $coverage_cutoffs;
system("cd $vfwd; $vfwd/generate_snp_properties.pl --manifest $manifest");
if ($?) \{
  exit 1;
\}

