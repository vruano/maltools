package MalariaGEN::AGV::Engines::Scratch::Sanger;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);
use File::stat;
use File::Path;
use File::Basename qw(dirname);

use Moose;

extends 'MalariaGEN::AGV::Engines::Scratch::Directory';

has 'low_space_threshold' => (is => 'rw', isa => 'Int', default => 100 * 1000 * 1000);

sub _irsync {
  my ($self,$source,$dest) = @_;
  File::Path::make_path(dirname($dest));
  -d dirname($dest) or die "cannot create destination directory " . dirname($dest) . "\n";
  `irsync -s i:$source $dest`;
  $? and die "could not synchronize irods '$source' with local file '$dest' $? $! $@";
}

sub _isize {
  my ($self,$df) = @_;
  my $result = `ils -l $df | awk '{print \$4}' | head -n 1`; chomp $result;
  $? and die "could not resolve file size for irods file $df, does it actually exists?";
  return $result;
}

sub _iget {
  my ($self,$df,$sf) = @_;
  my $df_size = $self->_isize($df);
  unless (-e $sf) {
    `iget $df $sf`;
    $? and die "did not manage to import $df from irods error code $? $! $@";
  }
  else {
    defined($df_size) or die "could not resolve file size for irods file $df, does it actually exists?";
    my $sf_stats = stat($sf);
    if ($sf_stats->size() != $df_size) {
      `iget -f $df $sf`; 
      $? and die "did not manage to re-import $df from irods error code $? $! $@";
    }
    else {
      print STDERR "Reusing exported irods file in '$sf' with the same size as in the repository\n"
    }
  }
  if (stat($sf)->size() != $df_size) {
    die "mismatch between irods size file $df_size and local copy " . stat($sf)->size();
  }
}

override '_import_type' => sub {
  my ($self,$file,$type,%args) = @_;
  $file = catfile(@$file) if ref($file) eq "ARRAY";
  my $df = $self->_data_file($file);
  my $sf = $self->_scratch_file($file);
  if ($df =~ /^\/seq/) { #irods!!!
    $self->_irsync($df,$sf);
    return $sf;
  }
  -e $df || $args{make_path} or die "data file to be imported does not exists '$df'";
  unless (-e $df) {
    $self->make_path($sf);
    return $sf;
  }
  my $stat_sf = stat($sf);
  my $stat_df = stat($df); 
  $self->_compat_filetype($stat_sf,$stat_df,$type) 
     || die "either source or destination are not compatible"
      if defined($stat_sf);
  if (($args{force} || 0) || !defined($stat_sf) 
     || $stat_sf->mtime < $stat_df->mtime) { 
   my $sf_dir = dirname($sf);
   File::Path::make_path($sf_dir);
   $self->_rcopy($df,$sf);
  }
  return $sf;
};

override 'quota' => sub {
  my $self = shift;
  my $root = $self->root;
  my $lfs_quota = `lfs quota $root`;
  my @lines = split(/\n/,$lfs_quota);
  my $line = $lines[3] or die "could not resolve the lfs quota for $root";
  chomp $line;
  my @fields = split(/\s+/,$line);
  @fields = grep { defined($_) && length("" . $_) > 0 } @fields;
  # Needed because some times values are give within [] brakets indicating inaccuracy:
  @fields = map { ($_ && $_ =~ /^\[(.*)\]$/) ? $1 : $_ } @fields;
  my ($kb_used,$kb_quota,$kb_limit,$kb_grace,$f_used,$f_quota,$f_limit,$f_grace) = @fields;
  my %result = (
     kb_used => $kb_used,
     kb_limit => $kb_limit,
     kb_remaining => $kb_limit - $kb_used,
  ); 
  my $low_space_thr = $self->low_space_threshold();
  if ($result{kb_remaining} < $low_space_thr) {
    print STDERR sprintf("WARNING!!! less than %.2f GB of space left in $root\n",$result{kb_remaining}/1000000); 
  }
  return wantarray ? %result : \%result;
};

1;

__END__
