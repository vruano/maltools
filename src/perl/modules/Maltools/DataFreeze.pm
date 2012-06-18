package Maltools::DataFreeze;

use JSON::XS;

use Moose;
use IO::File;
use IO::Pipe;

has 'file' => (is => 'ro', isa => 'Maybe[Str]', default => undef);
has 'fh' => (is => 'ro', isa => 'FileHandle', lazy => 1, builder => '_build_fh');
has '_already_read' => (is => 'rw', isa => 'Bool', default => 0);
has '_content' => (is => 'ro', isa => 'Any', lazy => 1, init_arg => 'content',
                   writer => '_set_content' , default => sub { {} });

sub BUILD {
  my $self = shift;

  unless (defined $self->file || defined $self->fh || defined $self->_content) {
    $self->_set_content({ samples => {}});
  } 
  $self->_read_freeze_file;
}

sub _build_fh {
  my $self = shift;
  my $file = $self->file;
  my $result;
  if ($file =~ /gz$/) {
    $result = IO::Pipe->new();
    $result->reader("gunzip -c < $file") or die "could not open the compressed freeze file '$file'";
  }
  else {
    $result = IO::File->new($file,'r') or die "could not open the compressed freeze file '$file'";
  }
  return $result;
}

sub _read_freeze_file {
  my $self = shift;
  return if ($self->_already_read);
  $self->_already_read(1);
  my $fh = $self->fh;
  my $json = JSON::XS->new();
  while (my $line = <$fh>) { $json->incr_parse($line) };
  $fh->close;
  $self->_set_content($json->incr_parse); 
}

sub find_sample {
  my ($self,$name) = @_;

  my $content = $self->_content();
  my $samples = $content->{samples};
  my $target = $samples->{$name};
  return $target;
}

sub find_sample_lane_and_files {
  my ($self,$name) = @_;
  my $sample_freeze = $self->find_sample($name);
  return undef unless $sample_freeze;
  my $lanes = $sample_freeze->{lanes};

  my %result = ();
  foreach my $lane (@$lanes) {
     my $lane_name = $lane->{run} ."." . $lane->{lane};
     my $files = $lane->{files} || [];
     @$files = grep { $_->{species_filter} eq "nonhuman" } @$files;
     next unless $#$files >= 0;
     die "multiple files available for the same lane '$lane_name' on sample '$name' " if $#$files >= 1;
     @$files = map { $_->{path} .  $_->{name} } @$files;
     $result{$lane_name} = $$files[0]; 
  }
  return wantarray ? %result : \%result;
}

1;
