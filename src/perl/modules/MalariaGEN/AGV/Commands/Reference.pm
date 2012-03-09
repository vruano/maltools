package MalariaGEN::AGV::Commands::Reference;

use Moose;
extends 'MalariaGEN::AGV::Command';

use POSIX;
use MalariaGEN::AGV::Config qw(reference_config);
use Getopt::Long qw(:config pass_through);

use IO::File;

has 'rc' => (is => 'ro', lazy => 1 , lazy => 1, default => sub { reference_config() });
has 'reference_file' => (is =>  'rw', isa => 'Str', lazy => 1, builder => '_build_reference_file' );
has 'reference_index_file' => (is => 'rw', isa => 'Str', lazy => 1, builder => '_build_reference_index_file' );
has 'interval_list_file' => (is => 'rw', isa => 'Maybe[Str]', lazy => 1, builder => '_build_interval_list_file' );
has 'sequences' => (is => 'rw', isa => 'Maybe[ArrayRef[Str]]', default => sub { undef });

sub help_summary {
   return "perform operation and queries upon the reference";
}

sub help_text {
   return "reference:\n\t" .
          $_[0]->help_summary ."\n".
          "Syntaxis:\n" .
          $_[0]->cl_name . " reference [query-name|command] [args]\n" .
          "Common arguments:\n" .
          "\t-r <file> : the name of the reference .fa file, we assume that the index is found in .fa.fai\n" .
          "\t\tIf none is provided it will take from the configuration\n" .
          "\t-i <file> : the index file\n" .
          "\t\tIf none is provided it will append .fai to the provided reference (-r) or just take it from the configuration\n" .
          "Queries:\n" .
          "\tintervals: generates a list of region interval suitable for dividing work in locus independent analyses\n" .
          "\t\tagv reference intervals [-n target-number-of-intervals| -s region-target-size]\n";
}

sub execute {
  my ($self,$samtrak_site,$params) = @_;
  local @ARGV = @{$params->{arguments}};
  my $reference = undef;
  my $reference_index = undef;
  GetOptions("reference|r=s" => \$reference, "index|i=s" => \$reference_index);
  $reference ||= $self->rc->file_name();
  $reference_index ||= $reference . ".fai";
  -f $reference or return $self->error_return("there reference file '$reference' is not reachable or not a regular file");
  -f $reference_index or return $self->error_return("there reference index '$reference_index' is not reachable or not a regular file");
  
  $self->reference_file($reference);
  $self->reference_index_file($reference_index);
  my $cmd = shift @ARGV;
  return $self->error_return("you must provide a command") unless defined $cmd;
  if ($cmd eq 'intervals') {
    return $self->intervals(@ARGV);
  }
  else {
    return $self->error_return("unknown command '$cmd'");
  }
}

sub intervals {
  my $self = shift;
  local @ARGV = @_;
  my $number = undef;
  my $size = undef;
  my $with_names;
  GetOptions('number|n=i' => \$number, 'size|s=i' => \$size, 'with-names|wn!' => \$with_names);
  my $intervals;
  unless (defined $number || defined $size) {
    if (defined $self->interval_list_file) {
      $intervals = $self->default_intervals();
    }
    else {
      return $self->error_return("you must provide at least one of the following: number of intervals (-n) or targeted size (-s)") unless defined $number || defined $size;
    }
  }
  else {
    $size = 1000000000 if !$size || $size <= 0;
    $number = 1 if !$number || $number <= 0;
    $intervals = $self->calculate_intervals($number,$size) or return $self->error_return("could not recover construct intervals to the given specifications");
  }
  $self->print_out_intervals($intervals,$with_names);
  return $self->ok_return();
}

sub default_intervals {
  my ($self) = @_;
  
  my $fh = IO::File->new($self->interval_list_file,'r') or return undef;
  my @result = ();

  while (my $line = <$fh>) {
    next if $line =~ /^(@|#)/ || $line !~ /\S/;
    chomp $line;
    my ($chr,$start,$stop,$sense,$name) = split(/\t/,$line);
    push @result, [ $chr,$start,$stop,$name];    
  }
  if ($self->sequences) {
     my %sequences = map { $_ => 1 } @{$self->sequences};
     @result = grep { exists ($sequences{$$_[0]}) } @result ;
  }

  return wantarray ? @result : \@result;
}

sub read_index {
  my $self = shift;
  my @result = ();
  my $fh = IO::File->new($self->reference_index_file,'r') or return undef;
  while (my $line = <$fh>) {
     chomp $line;
     push @result, [ split (/\t/,$line) ];
  }
  $fh->close;
  return \@result;
}


sub calculate_intervals {
  my $self  = shift;
  my $number = shift;
  my $size = shift;
  my $seqs = $self->read_index() or return $self->error_return("could not recover any sequence from index file '" . $self->reference_index_file . "'");
  if ($self->sequences) {
     my %sequences = map { $_ => 1 } @{$self->sequences};
     $seqs = [ grep { exists ($sequences{$$_[0]}) } @$seqs ];
  } 

  my $genome_length = 0;
  foreach my $s (@$seqs) { $genome_length += $$s[1] };
 
  my $size_by_number = ceil($genome_length / $number);
  $size = $size_by_number if ($size > $size_by_number);
  my $number_by_size = ceil($genome_length / $size); 
  $number = $number_by_size if ($number < $number_by_size);
  my @result = ();
  my $index = 0;

  foreach my $s (@$seqs) {
    my $sequence_length = $$s[1];
    my $start = 1;
    while ($start < $sequence_length) {
      my $end = $start + $size;
      $end = $sequence_length if $end > $sequence_length;
      push @result, [ $$s[0], $start, $end, "interval_" . $index++ ];
      $start = $end + 1; 
    }
  }
  return \@result;
}

sub print_out_intervals {
   my ($self,$intervals,$with_names) = @_;
   foreach my $inv (@$intervals) {
     print STDOUT $$inv[0] . ":" . $$inv[1] . "-" . $$inv[2] . ($with_names ? "\t$$inv[3]" : "") . "\n";
   }
}

sub _build_reference_file {
  my $self = shift;
  my $rc = $self->rc();
  return $rc->file_name();
}

sub _build_reference_index_file {
  my $self = shift;
  my $reference_file = $self->reference_file;
  $reference_file .= ".fai";
  return -e $reference_file ? $reference_file : undef;
}

sub _build_interval_list_file {
   my $self = shift;
   my $reference_file = $self->reference_file;
   $reference_file =~ s/((\.fa)|(\.fasta))$/.interval_list/;
   return -f $reference_file ? $reference_file : undef;
}

1;
