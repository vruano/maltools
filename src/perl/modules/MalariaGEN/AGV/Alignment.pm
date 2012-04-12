package MalariaGEN::AGV::Alignment;

use Moose;

use File::Temp qw(tempdir);
use Bio::DB::Sam;
use MalariaGEN::AGV::Reference;
use MalariaGEN::AGV::Alignment::ReadGroup;
use File::Spec::Functions qw(catfile);

has file => (is => 'ro', isa => 'Str', required => 1);
has reference => (is => 'ro', isa => 'Maybe[MalariaGEN::AGV::Reference]');
has 'header_text' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_header_text');
has 'header_lines' => (is => 'ro', isa => 'ArrayRef[Str]', lazy => 1, builder => '_build_header_lines');
has 'dictionary' => (is => 'ro', isa => 'ArrayRef', lazy => 1, builder => '_build_dictionary');
has '_db' => (is => 'ro', isa => 'Any', lazy => 1, builder => '_build_db');
has read_groups => (is => 'ro', isa => 'ArrayRef[MalariaGEN::AGV::Alignment::ReadGroup]', lazy => 1, builder => '_build_read_groups');
has _fastq_outdir => (is => 'rw', isa => 'Str', init_arg => undef, lazy => 1, builder => '_build_fastq_outdir');

sub order {
   my $self = shift;
   my @header_lines = @{$self->header_lines};
   return "unsorted" unless $#header_lines >= 0;
   return lc($1) if $header_lines[0] =~ /\tSO:(\S+)/;
   return "unsorted";
}

sub _build_fastq_outdir {
  my $self = shift;
  my $result = tempdir( "FastqXXXXX", CLEANUP => 0 );
  mkdir $result;
  $self->to_fastq(outdir => $result);
  return $result;
}

sub _forward_fastq_file {
  my $self = shift;
  my $group = shift;
  my $id = $group->id;
  my $outdir = $self->_fastq_outdir;
  return catfile($outdir,$id . "_1.fastq") if -e catfile($outdir,$id . "_1.fastq");
  return catfile($outdir, $group->platform_unit . "_1.fastq") if $group->platform_unit;
  
}

sub _reverse_fastq_file {
  my $self = shift;
  my $group = shift;
  my $id = $group->id;
  my $outdir = $self->_fastq_outdir;
  return catfile($outdir,$id . "_2.fastq") if -e catfile($outdir,$id . "_2.fastq");
  return catfile($outdir, $group->platform_unit . "_2.fastq") if $group->platform_unit;
  
}

sub to_fastq {
   my $self = shift;
   my %args = @_;
   if ($args{outdir}) {
     my $file = $self->file;
     my $outdir = $args{outdir};
     `mkdir -p $outdir`;
     die "could not create output directory $outdir (code=$?)" if $?;
     `picard RevertSam COMPRESSION_LEVEL=0 INPUT=$file OUTPUT=/dev/stdout | picard SamToFastq INPUT=/dev/stdin OUTPUT_DIR=$outdir OUTPUT_PER_RG=true`;
     die "error generating fastq (code=$?)" if $?;
     foreach my $rg (@{$self->read_groups}) {
       my $pu = $rg->platform_unit;
       my $id = $rg->id;
       next unless defined $pu;
       my $f1 = $pu . '_1.fastq';
       my $f2 = $pu . '_2.fastq';
       my $new_f1 = $f1;
       $new_f1 =~ s/$pu/$id/;
       my $new_f2 = $f2;
       $new_f2 =~ s/$pu/$id/;
       $f1 =~ s/#/_/g;
       $f2 =~ s/#/_/g;
       rename catfile($outdir,$f1), catfile($outdir,$new_f1) if -e catfile($outdir,$f1);
       rename catfile($outdir,$f2), catfile($outdir,$new_f2) if -e catfile($outdir,$f2);
     }
   }
   elsif ($args{forward} || $args{reverse}) {
     my $forward = $args{forward} || '/dev/null';
     my $reverse = $args{reverse} || '/dev/null';
     my $file = $self->file;
     `picard SamToFastq INPUT=$file FASTQ=$forward SECOND_END_FASTQ=$reverse`;
     die "error generating fastq (code=$?)" if $?;
   }
   else {
     die "you need to indicate an output directory or a forward/reverse file name";
   }
}

sub valid_algorithm {
   my $self = shift;
   $#_ >= 0 or die "you need to provide an algorithm to validate against";
   my $alg;
   if ($#_ == 0) {
     $alg = ref($_[0]) eq "" ? {name => shift} : shift;
   }
   else {
     $alg = { @_ };
   }
   
   my $alg_name = $$alg{name} or die "no algorithm name was provided";
   
   my $header_lines = $self->header_lines;
   my @pg_lines = grep { index($_,'@PG') == 0 } @$header_lines;
   return 0 unless $#pg_lines >= 0;
   @pg_lines = grep { $_ =~ /\tPN:\Q${alg_name}\E(\t.*)?$/ } @pg_lines; 
   return 0 unless $#pg_lines >= 0;
   return 1 unless exists $$alg{version};
   my $version = $$alg{version};
   @pg_lines = grep { $_ =~ /\tVN:\Q${version}\E(\t.*)?$/ } @pg_lines;
   return scalar(@pg_lines);
}

sub is_sorted {
   return $_[0]->order ne "unsorted";
}

sub valid_reference {
   my $self = shift;
   my $ref = shift;
   my $my_dictionary = $self->dictionary;
   return $ref->same_dictionary($my_dictionary);
}

sub sorted_by_coordinates {
   return $_[0]->order =~ /^coord/;
}

sub _build_db {
   my $self = shift;
   my $file = $self->file or die "no file was provided";
   return Bio::DB::Sam->new(-bam => $file);
}

sub _build_header_text {
   my $self = shift;
   my $db = $self->_db;
   return $db->header->text;
}

sub _build_read_groups {
  my $self = shift;
  my $hl = $self->header_lines;
  return [ map { MalariaGEN::AGV::Alignment::ReadGroup->new(line => $_, alignment => $self) } grep { index($_,'@RG') == 0 } @$hl ];
}

sub _build_header_lines {
  my $header_text = $_[0]->header_text;
  my @lines = split(/\n/m,$header_text);
  for (my $i = 0; $i <= $#lines; $i++) {
    chomp $lines[$i];
  }
  return \@lines;
}

sub _build_dictionary {
  my @header_lines = @{$_[0]->header_lines};

  my @result = ();

  foreach my $hl (@header_lines) {
    next unless $hl =~ /^\@SQ/;
    my @c = split(/\t/,$hl);
    my %c = (map { (split(/:/,$_,2))[0..1]} @c);
    push @result, { name => $c{SN}, length => $c{LN}, md5 => $c{M5} };
  }
  
  if (scalar(@result) == 0 && defined $_[0]->reference) {
     return $_[0]->reference->dictionary;
  }
  return \@result;
}




1;
