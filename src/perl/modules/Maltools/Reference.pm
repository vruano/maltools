package Maltools::Reference;

use Moose;
use POSIX;
use IO::File;

with 'MooseX::Cloneable';

use overload
   fallback => 1,
   '""' => sub { $_[0]->sequence_file };

has 'sequence_file' => (is => 'ro', isa => 'Str', required => 1);
has 'index_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_index_file_builder');
has 'dictionary_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_dictionary_file_builder');
has 'annotation_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_annotation_file_builder');
has 'chromosomes' => (is => 'ro', isa => 'ArrayRef[Str]', lazy => 1, builder => '_build_chromosomes');
has 'dictionary' => (is => 'ro', isa => 'ArrayRef', lazy => 1, builder => '_build_dictionary');
has 'interval_file' => (is => 'ro', isa => 'Maybe[Str]', lazy => 1, builder => '_build_interval_file');
has 'possible_snps_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_possible_snps_file');
has 'probable_snps_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_probable_snps_file');
has 'uniqueness_file' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_uniqueness_file');


around BUILDARGS => sub {
   my $orig = shift;
   my $class = shift;

   if (@_ == 1 && !ref($_[0])) {
      return $class->$orig(sequence_file => shift)
   }
   my %args;
   if (@_ == 1 && ref($_[0]) eq "HASH") {
      %args = %{$_[0]};
   }
   elsif (@_ == 1 && ref($_[0]) eq "ARRAY") {
      %args = @{$_[0]};
   }
   else {
      %args = @_;
   }

   $args{sequence_file} ||= $args{sequences} || $args{fasta} || $args{fa} if (exists($args{sequences}) || exists($args{fa}));
   $args{index_file} ||= $args{index} || $args{fai} if (exists($args{index}) || exists($args{fai}));
   $args{dictionary_file} ||= $args{dictionary} if (exists($args{dictionary}) && ref($args{dictionary}) eq "");
   $args{dictionary_file} ||= $args{dict} if (exists($args{dict}) && ref($args{dict}) eq "");
   $args{dictionary} ||= $args{dict} if (exists($args{dict}) && ref($args{dict}) eq "ARRAY");
   $args{annotation_file} ||= $args{annotation} || $args{gff} if (exists($args{annotation}) || exists($args{gff}));

   delete @args{qw(sequence fasta index fa fai dict annotation gff)};
   return $class->$orig(%args);
};

sub exists {
  my $self = shift;
  return -e $self->sequence_file;
}

sub location {
  my $self = shift;
  my $chr = shift;
  my $pos = shift || 1;
  require Maltools::Reference::Location;
  return Maltools::Reference::Location->new(reference => $self, sequence => $chr, position => $pos);
}

sub _annotation_file_builder {
   my $self = shift;
   my $result = $self->sequence_file;
   $result =~ s/\.fa(sta)?(\.gz)?//;
   $result .= ".gff";
   return $result if -f $result;
   return undef;
}

sub _index_file_builder {
   my $self = shift;
   my $result = $self->sequence_file . ".fai";
   die "unreachable reference index file '$result' or not a file" unless -f $result;
   return $result if -f $result;
   return undef;
}

sub _dictionary_file_builder {
   my $self = shift;
   my $ref_file = $self->sequence_file;
   $ref_file =~ s/\.fa(sta)?(\.gz)?$//;
   return $ref_file . ".dict";
}

sub _build_uniqueness_file {
   my $self = shift;
   my $ref_file = $self->sequence_file;
   $ref_file =~ s/\.fa(sta)?(\.gz)?$//;
   return $ref_file . ".uq";
}


sub _build_possible_snps_file {
   my $self = shift;
   my $result = $self->sequence_file;
   $result =~ s/\.fa(sta)?(\.gz)?$/.snps/;
   return $result if -f $result;
   return $result . ".gz" if -f $result . ".gz";
   return undef;
}

sub _build_probable_snps_file {
   my $self = shift;
   my $result = $self->sequence_file;
   $result =~ s/\.fa(sta)?(\.gz)?$/.snps-mask.bed/;
   return $result if -f $result;
   return $result . ".gz" if -f $result . ".gz";
}

sub _build_interval_file {
   my $self = shift;
   my $result = $self->sequence_file;
   $result =~ s/\.fa(sta)?(\.gz)?$/.intervals/;
   return $result if -f $result;
   return $result . ".gz" if -f $result . ".gz";
   return undef;
}

sub _interval_file_contents {
   my $self = shift;
   my $ifc = $self->interval_file;
   return undef unless $ifc;
   die "interval file '$ifc' does not exists, is not a regular file, cannot be accessed or read" unless -r $ifc && -f $ifc;
   open my $ifc_fh , ($ifc =~ /.gz$/) ? "gunzip -c $ifc |" : $ifc
      or die "could not open or decompress interval file '$ifc' ";
   my @result = ();
   while (my $line = <$ifc_fh>) {
      next unless $line =~ /\S/;
      next if $line =~ /^\@/;
      chomp $line;
      push @result, [split (/\t/,$line)];
   }
   return wantarray ? @result : \@result;
}

sub _build_dictionary {
   my $self = shift;
   my $dict_file = $self->dictionary_file;
   -f $dict_file or die "dictionary file '$dict_file' is not accessible";
   my $fh = IO::File->new($dict_file,'r') or die "could not open dictionary file '$dict_file'";
   
   my @result = ();
   while (my $line = <$fh>) {
     chomp $line;
     my @c = split(/\t/,$line);
     next unless $#c >= 0;
     next unless $c[0] eq '@SQ';
     my %c = map { split(/:/,$_,2) } @c[1..$#c];
     push @result, {name => $c{SN}, md5 => $c{M5}, length => int($c{LN}) };
   }
   $fh->close();
   return \@result;
}

sub same_dictionary {
   my $self = shift;
   my $other = shift;
   my $mine = $self->dictionary;
   return 0 if ($#$other != $#$mine);
   return 0 unless $#$other >= 0;
   for (my $i = 0; $i <= $#$other; $i++) {
      my $o = $$other[$i];
      my $m = $$mine[$i];
      return 0 unless $o->{name} eq $m->{name};
      return 0 unless $o->{length} == $m->{length};
      return 0 unless defined($o->{md5}) == defined($m->{md5});
      return 0 unless $o->{md5} eq $m->{md5};
   }
   return 1;
}



sub _index_contents {
  my $self = shift;
  my $index_file = $self->index_file;
  my $index_fh = IO::File->new($index_file);
  my @result = ();
  while (my $line = <$index_fh>) {
    chomp $line;
    push @result, [split(/\t/,$line)];
  }
  return wantarray ? @result : \@result;
}

sub create_intervals {
  my $self = shift;
  my %args = @_;
  my $min_number = $args{min_number} || 1;
  my $max_size = $args{max_size} || 10000000000;
  my @base_intervals = $self->_create_base_intervals(%args);
  my ($total_length,$max_length) = $self->_create_intervals_total_length(@base_intervals);
  # calculate the actual number of targeted intervals based on size 
  my $size_by_number = ceil($total_length / $min_number);
  my $number_by_size = ceil($total_length / $max_size);
  $max_size = $size_by_number if ($max_size > $size_by_number);
  $min_number = ($min_number > $number_by_size) ? $min_number : $number_by_size;
  # just return the base collection of intervals if this satisfies the requirements:
  if ( (! $args{max_size} || ($args{max_size} >= $max_length) ) 
       && scalar(@base_intervals) >= $min_number) {
    return wantarray ? @base_intervals : \@base_intervals;
  }

  # otherwise do the job: 
  my @result = ();
  foreach my $itv (@base_intervals) {
    my $start = $$itv[1];
    while ($start <= $$itv[2]) {
       my $end = $start + $max_size - 1;
       $end = $$itv[2] if $$itv[2] < $end;
       push @result, [ $$itv[0], $start, $end ];
       $start = $end + 1;
    }
  }
  return wantarray ? @result : \@result;
}

sub _create_intervals_total_length {
  my $self = shift;
  my $total_length = 0;
  my $max_length = 0;
  foreach my $itv (@_) {
    my $len = $$itv[2] - $$itv[1] + 1;
    $max_length = $len if ($len >= $max_length);
    $total_length += $len;
  }
  return ($total_length,$max_length);
}

sub has_interval_file {
   return defined $_[0]->interval_file;
}

sub _create_base_intervals {
  my $self = shift;
  my %args = @_;
  my @intervals = @{$args{intervals} ||[]};
  if (scalar(@intervals) == 0) {
     if ($self->has_interval_file) {
       @intervals = map { [@$_[0..2]] }  $self->_interval_file_contents();
     }
     else {
       @intervals = map { [$$_[0],1,$$_[1]]  } $self->_index_contents();
     }
  } 
  for (my $i = 0; $i <= $#intervals; $i++) {
     my $itv = $intervals[$i];
     if (ref($itv) eq "") {
        $itv =~ /(\S+):(\d+)\-(\d+)/ or die "incorrect internal spec format '$itv'";
        die "incorrect interval range start ($2) > end ($3) in '$itv' " if ($2 > $3);
	$intervals[$i] = [ $1, $2 , $3 ];
     }
  }
  return @intervals;
}

sub calculate_intervals {
  my $self  = shift;
  my $number = shift;
  my $size = shift;
  my $seqs = shift;

  my $genome_length = 0;
  foreach my $s (@$seqs) { $genome_length += $$s[1] };

  my $size_by_number = ceil($genome_length / $number);
  $size = $size_by_number if ($size > $size_by_number);
  my $number_by_size = ceil($genome_length / $size);
  $number = $number_by_size if ($number < $number_by_size);
  my @result = ();

  foreach my $s (@$seqs) {
    my $sequence_length = $$s[1];
    my $start = 1;
    while ($start < $sequence_length) {
      my $end = $start + $size;
      $end = $sequence_length if $end > $sequence_length;
      push @result, [ $$s[0], $start, $end ];
      $start = $end + 1;
    }
  }
  return \@result;
}

sub write_intervals_list_file {
  my ($self,%args) = @_;
  my $fh = $args{fh};
  unless ($fh) {
    if ($args{file}) {
      $fh = IO::File->new($args{file},'w') or die "could not open intervals list file '" . $args{file} ."' for writing";
    }
    else {
      die "no file handle nor file name provided";
    }
  }
  print $fh join("\t",'@HD','VN:1.0','SO:coordinate'),"\n";
  foreach my $seq ($self->_index_contents) {
    print $fh join("\t",'@SQ','SN:' . $$seq[0],'LN:' . $$seq[1]),"\n";
  }
  my $next_index = 0;
  foreach my $interval ($self->create_intervals(%args)) {
    print $fh join("\t",@$interval,'+',"interval_" . $next_index++),"\n"; 
  }
  $fh->close unless $args{fh};
}

sub _build_chromosomes {
  my ($self) = shift;

  my $rif = $self->index_file;
  
  if ($rif && -f $rif) {
    return $self->_build_chromosomes_from_index_file($rif);
  }

  my $rf = $self->sequence_file;

  if ($rf && -f $rf) {
    return $self->_build_chromosomes_from_sequence_file($rf);
  }

  die "there is no reference file specified";
}

sub _build_chromosomes_from_index_file {
  my ($self,$rif) = @_;

  my $fh = IO::File->new($rif,'r') or return [];
  my @result = ();
  while (my $line = <$fh>) {
    chomp $line;
    my ($chr) = split (/\t/,$line);
    push @result, $chr;
  }
  $fh->close();
  return \@result;

}

sub _build_chromosomes_from_sequence_file {
  my ($self,$rf) = @_;
  my $fh = IO::File->new($rf,'r') or return [];
  my @result = ();
  while (my $line = <$fh>) {
    chomp $line;
    next unless $line =~ /^>(\S+)/;
    push @result, $1;
  }
  $fh->close();
  return \@result;
}



1;
