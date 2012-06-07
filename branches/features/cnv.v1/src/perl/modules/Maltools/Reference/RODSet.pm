package Maltools::Reference::RODSet;

use Moose;
use Maltools::Reference::ROD;
use Cwd qw(getcwd);
use File::Spec::Functions qw(file_name_is_absolute rel2abs catfile);

has '_rods_by_name' => (
    is => 'ro',
    isa => 'HashRef[Maltools::Reference::ROD]',
    builder => '_build_rods_by_name',
    lazy => 1,
    traits => ['Hash'],
    handles => {
        get => 'get',
        names => 'keys',
        _add => 'set'
    },
);

has 'list_file' => (is => 'ro', isa => 'Str', required => 0);
has 'search_path' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_search_path');

sub _build_search_path {
  my $self = shift;
  return getcwd();
}


sub _build_rods_by_name {
  my $self = shift;
  my $list_file = $self->list_file;
  return {} unless ($list_file);
  -e $list_file or die "cannot reach the ROD list file '$list_file'";
  -f $list_file or die "the ROD list file '$list_file' is in fact not a ROD file";
  open my $fh , $list_file or die "could not open '$list_file' read-only";
  my %result = ();
  while (my $line = <$fh>) {
     next unless $line =~ /\S+/;
     chomp $line;
     my ($name,$type,$file) = split(/\t/,$line);
     if (!file_name_is_absolute($file)) {
        $file = rel2abs($file,$self->search_path);
     }
     $result{$name} = Maltools::Reference::ROD->new(type => $type, file => $file);
  }
  return \%result;
}


1;