package Maltools::Manifest::Pipeline::Sample;

use Moose;

use Cwd qw(realpath);
require Maltools::Manifest::Pipeline::DataSource;
require Maltools::Manifest::Pipeline;

has manifest => (is => 'ro', isa => 'Maltools::Manifest::Pipeline', required => 1);
has oxfordCode => (is => 'rw', isa => 'Str');
has samtrakId => (is => 'rw', isa => 'Str');
has dataSources => (is => 'ro', isa => 'ArrayRef[Maltools::Manifest::Pipeline::DataSource]', default => sub { [] });


around BUILDARGS => sub {
  my $orig = shift;
  my $class = shift;
  my %args = @_;
  $args{dataSources} = [ map { blessed($_) ? $_ : Maltools::Manifest::Pipeline::DataSource->new(%$_) } @{$args{dataSources} || []} ];
  return $class->$orig(%args);
};

sub has_single_data_source {
    return scalar(@{$_[0]->dataSources}) == 1;
}

sub has_inout_data_source {
    my $self = shift;
    return 0 unless $self->has_single_data_source;
    my $data = ${$self->dataSources}[0];
    my $url = $data->url;
    $url = $1 if $url =~ /^file:\/\/?(\S+)$/;
    return 0 if ($url =~ /^[a-zA-Z\-]+:/);
    
    my $real_ds = realpath($url);
    my $real_ap = realpath($self->manifest->get_path('alignments'));

    return index($real_ds,$real_ap) == 0;
}



1;