package MalariaGEN::AGV::Manifest::Pipeline::DataSource;

use Moose;

use File::Basename qw(basename);

has 'url' => (is => 'ro', isa => 'Str', required => 1);
has 'md5' => (is => 'rw', isa => 'Str', predicate => 'has_md5');

sub file_basename {
    my $self = shift;
    my $url = $self->url;
    $url =~ s/^([^:]+:\/*)//;
    return basename($url);
}

sub full_filename {
    my $self = shift;
    my $url = $self->url;
    $url =~ s/^[^:]+:\/*//;
    $url = "/" . $url;
    return $url; 
}

1;