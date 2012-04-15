package MalariaGEN::AGV::Engines::Scratch::Array;

use strict;
use warnings;
use Cwd qw(realpath);
use Carp;
use Moose;

with 'MalariaGEN::AGV::Engines::Scratch';


has _last_best_time => (is => 'rw', default => 0);
has _last_best_area => (is => 'rw', default => undef);
has areas => (is => 'ro', isa => 'ArrayRef[MalariaGEN::AGV::Engines::Scratch]', default => sub { [] });
has exclusions => (is => 'ro', isa => 'ArrayRef[Str]', default => sub { [] });

# Indicate what is the best scratch are to use for a given set of parameters:
# 1. checks whether the file or directory involved directly points to an area.
# 2. checks the argument close_to just the same as above.
# 3. checks whether the file seems to exist already in some area. 
# 4. finally checks the quota remaining in each area and returns the emptiest one.
sub _best_area {
  my $self = shift;
  my %args = @_;
  my @areas = @{$self->areas()};
  my $file_or_dir = $args{file} || $args{dir} || $args{file_or_dir};
  # If fi
  if (defined $file_or_dir) {
    my @file_areas = grep { index($file_or_dir,$_->root) >= 0 } @areas;
    return $file_areas[0] if $#file_areas >= 0;
  }
  if ($args{close_to}) {
    my $close_to = $args{close_to};
    my @close_areas = grep { index($close_to,$_->root) >= 0 } @areas;
    return $close_areas[0] if $#close_areas >= 0;
  }
  if (defined $file_or_dir) {
    my @file_areas = grep { $_->in_scratch($file_or_dir) } @areas;
    return $file_areas[0] if $#file_areas >= 0;
  }
  if (time() - $self->_last_best_time < 10) {
    return $self->_last_best_area;
  } 
  my @quotas = map { scalar($_->quota) } @areas;
  my @indexes =  sort { $quotas[$b]->{kb_remaining} <=> $quotas[$a]->{kb_remaining} } (0..$#areas);
  my $result = $areas[$indexes[0]];
  $self->_last_best_time(time());
  $self->_last_best_area($result);
  return $result;
}

sub in_scratch {
  my $self = shift;
  my $target = shift;
  return grep { $_->in_scratch($target) } @{$self->areas()};
}

sub tempdir {
 my ($self,%args) = @_;
 my $best = $self->_best_area(%args);
 return $best->tempdir(%args);
}

sub quota {
 my ($self,%args) = @_;
 my $best = $self->_best_area(%args);
 return $best->quota();
}

sub tempfile {
  my ($self,%args) = @_;
  my $best = $self->_best_area(%args);
  return $best->tempfile(%args);
}

sub to_scratch {
  my ($self,%args) = @_;
  my $file_or_dir = $args{file_or_dir} || $args{file} || $args{dir};
  defined $file_or_dir or confess "must provide a file or directory name";
  -e $file_or_dir || $args{make_path} or confess "must provide an existing file or directory name '$file_or_dir'";
  if ($self->is_excluded(%args)) {
    return realpath($args{file_or_dir} || $args{file} || $args{dir});
  }
  my $best = $self->_best_area(%args);
  return $best->to_scratch(%args);
}

sub is_excluded {
  my ($self,%args) = @_;
  my $file_or_dir = $args{file_or_dir} || $args{file} || $args{dir};
  $file_or_dir = realpath($file_or_dir) || $file_or_dir;
  foreach my $excl (@{$self->exclusions}) {
    return 1 if index($file_or_dir,realpath($excl)) != -1;
  } 
  return 0;
}

sub from_scratch {
  my ($self,%args) = @_;
  if ($self->is_excluded(%args)) {
    return realpath($args{file} || $args{dir});
  }
  my $best = $self->_best_area(%args);
  return $best->from_scratch(%args);
}


1;
