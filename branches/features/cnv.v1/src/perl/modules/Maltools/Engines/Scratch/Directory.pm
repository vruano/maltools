package Maltools::Engines::Scratch::Directory;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);
use File::stat;
use File::Rsync;
use File::Path;
use File::Basename qw(dirname);
use File::Temp;
use Cwd qw(realpath);
use Fcntl qw(:mode);

use Moose;

with 'Maltools::Engines::Scratch';

has 'temp_root' => ( is => 'ro' , isa => 'Str' , lazy => 1, default => sub {my $res = catfile($_[0]->root,"tmp"); 
              File::Path::make_path($res); return $res; });
has '_rsync' => (is => 'ro', lazy => 1, default => sub { File::Rsync->new({archive => 1, compress => 0}) });
has 'root' => (is => 'ro', isa => 'Str', required => 1, writer => '_set_root' );

sub BUILD {
  my $self = shift;
  my $root = $self->root;
  my $real_root = realpath($root);
  $self->_set_root($real_root) if $root ne $real_root;
}

sub tempdir {
  my ($self,%args) = @_;
  my $template = exists $args{pattern} ? $args{pattern} : "XXXXXXXX"; 
  return File::Temp::tempdir($template,DIR => catfile($self->temp_root));
}

sub tempfile {
  my ($self,%args) = @_;
  my $unlink = $args{unlink} || 0;
  delete $args{unlink};
  my $template = exists $args{template} ? $args{template} : "XXXXXXXX";
  my $dir = ($args{dir} || catfile($self->temp_root));
  -e $dir or File::Path::make_path($dir);
  my ($fh,$fn) =  File::Temp::tempfile($template, DIR => $dir);
  $fh->close;
  unlink($fn) if $unlink; 
  return $fn;
}

sub make_path {
  my ($self,$path) = @_;
  $path = catfile($self->root,$path);
  File::Path::make_path(dirname($path));
}

sub to_scratch {
  my ($self,%args) = @_;
  my $file = $args{file};
  if (defined($file)) {
    return $self->_import_file($file,%args);
  }
  my $dir = $args{dir};
  if (defined($dir)) {
    return $self->_import_dir($dir,%args);
  }
  my $fod = $args{file_or_dir};
  my $afod = $self->_data_file($fod);
  if (defined($fod)) {
    if (-f $afod) {
      return $self->_import_file($fod,%args);
    }
    elsif (-d $afod) {
      return $self->_import_dir($fod,%args);
    }
    elsif (! -e $afod) {
      return $self->_import_file($fod,%args) if ($args{make_path});
      die "cannot find file to import '$fod' ";
    }
    else {
      die "unsupported file type for import '$fod'";
    }
  }
}

sub from_scratch {
  my ($self,%args) = @_;
  my $file = $args{file};
  if (defined($file)) {
    return $self->_export_file($file,%args);
  }
  my $dir = $args{dir};
  if (defined($dir)) {
    return $self->_export_dir($dir,%args);
  }
  my $fod = $args{file_or_dir};
  my $afod = $self->_scratch_file($fod);
  if (defined($afod)) {
    if (-f $afod) {
      return $self->_export_file($fod,%args);
    }
    elsif (-d $afod) {
      return $self->_export_dir($fod,%args);
    }
    else {
      die "cannot export file due to wrong type '$fod' '$afod'";
    }
  }

}

sub _export_type {
  my ($self,$file,$type,%args) = @_;
  $file = catfile(@$file) if ref($file) eq "ARRAY";
  my $df = $self->_data_file($file);
  die "irod files cannot exported\n" if ($df =~ /^\/seq/);
  my $sf = $self->_scratch_file($file);
  -e $sf or die "data file to be imported does not exists '$df'";
  my $stat_sf = stat($sf);
  my $stat_df = stat($df); 
  $self->_compat_filetype($stat_sf,$stat_df,$type) 
     || die "either source or destination are not compatible"
      if defined($stat_df);
  if (($args{force} || 0) || !defined($stat_df) 
     || $stat_sf->mtime > $stat_df->mtime) { 
   my $df_dir = dirname($sf);
   File::Path::make_path($df_dir);
   $self->_rcopy($sf,$df);
  }
  return $df;
}

sub _compat_filetype {
  my ($self,$file1,$file2,$type) = @_;
  my $mode1 = defined($file1) ? $file1->mode : undef;
  my $mode2 = defined($file2) ? $file2->mode : undef;
 
  if (defined $type) {
    return  0 if defined $mode1 && ! $type->($mode1);
    return  0 if defined $mode2 && ! $type->($mode2);
  }
  return 1 unless defined $mode1 &&  defined $mode2;
  
  my ($ftype1,$ftype2) = (S_IFMT($mode1), S_IFMT($mode2));
  return $ftype1 & $ftype2;
}

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

sub _import_type {
  my ($self,$file,$type,%args) = @_;
  $file = catfile(@$file) if ref($file) eq "ARRAY";
  my $df = $self->_data_file($file);
  my $sf = $self->_scratch_file($file);
  my $exists_df = -e $df;
  $exists_df || $args{make_path} or die "data file to be imported does not exists '$df'";
  unless ($exists_df) {
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
}


sub _export_file {
  my ($self,$file,%args) = @_;
  return $self->_export_type($file,\&S_ISREG,%args);
}

sub _export_dir {
  my ($self,$dir,%args) = @_;
  return $self->_export_type($dir,\&S_ISDIR,%args);
}


sub _import_file {
  my ($self,$file,%args) = @_;
  my $result = $self->_import_type($file,\&S_ISREG,%args);
  File::Path::make_path(dirname($result));
  return $result;
}

sub _import_dir {
  my ($self,$dir,%args) = @_;
  my $result = $self->_import_type($dir,\&S_ISDIR,%args);
  File::Path::make_path($result);
  return $result;
}

sub _data_file {
  my ($self,$file) = @_;
  $file = realpath($file) || $file;
  return $file;
}

sub _scratch_file {
  my ($self,$file) = @_;
  $file = realpath($file) || $file;
  return catfile($self->root,substr($file,1)) if ($file =~ /^\//);
  return catfile($self->root,$file);
}

sub _rcopy {
  my ($self,$source,$dest) = @_;
  $source .= '/' if (-d $source && $source !~ /\/$/);
  return if ($source eq $dest);
  $self->_rsync()->exec({ src => $source, dest => $dest}) or warn "rsync from '$source' to '$dest' failed $? $! $@\n";
}

sub quota {
  my $self = shift;
  my $root = $self->root;
  return 100000000;
  my @lines = `df -B 1024 $root`;
  $? and die "could not determine quota";
  
  my $last_line = $lines[$#lines];
  chomp $last_line;
  $last_line =~ /(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$/
    or die "error in df format '$last_line'";
  my ($total,$used,$free,$use_pc,$mount_point) = ($1,$2,$3,$4,$5);

  my %result = (kb_used => $used, kb_limit => $total, kb_remaining => $free);
  return wantarray ? %result : \%result;
}


sub in_scratch {
  my $self = shift;
  my $file = shift;
  my $sfile = $self->_scratch_file($file);
  return -e $sfile;
}

1;

__END__
