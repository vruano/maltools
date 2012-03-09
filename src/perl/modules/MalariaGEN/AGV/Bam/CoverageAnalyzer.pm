package MalariaGEN::AGV::Bam::CoverageAnalyzer;

use strict ;
use warnings ;
use JSON::XS;
use MalariaGEN::AGV::GffUtils::Iterator;
use Moose;
use POSIX qw(ceil floor sqrt abs);
use Scalar::Util qw(blessed);
use MalariaGEN::AGV::Config;

our @DEFAULT_PERCENTILES = qw(0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100);
our @DEFAULT_TYPES = qw(all exonic genic genic_non_coding intergenic non_coding featureless); 

has 'percentiles' => (is => 'ro', isa => 'ArrayRef', default => sub {[@DEFAULT_PERCENTILES]});
has 'gff' => (is => 'ro', required => 1);
has 'filters' => (is => 'ro', default => sub { [] });
has 'filter_all' => (is => 'ro', isa => 'Bool', default => 0);
has 'reference' => (is => 'ro', isa => 'MalariaGEN::AGV::Reference', default => sub { return MalariaGEN::AGV::Config->new()->get_reference(); } );

sub zero_sample_counter_group {
  my $self = shift;
  return { all => $self->zero_counter_group(),
           nuclear => $self->zero_counter_group(),
           map { $_ => $self->zero_counter_group() } @_,
  };
}

sub zero_counter_group {
 my $self = shift;
 return { map { $_ => $self->zero_counter() } @DEFAULT_TYPES};
}            

sub zero_counter {
 my $self = shift;
 { count => 0, sqtotal => 0, total => 0, max => 0, frequencies => [0],
   mean => undef, variance => undef, sd => undef,
   percentiles => undef, cumulative => [0] };
}

sub increase_counters {
  my ($self,$count_delta,$total_delta,$all_samples,$types,$counters) = @_;
  foreach my $cn (@$counters) {
    my $c = $$all_samples{$cn} || ($$all_samples{$cn} = $self->zero_counter_group());
    foreach my $t ( @$types ) {
      my $ct = $c->{$t};
      $ct->{count} += $count_delta;
      $ct->{total} += $total_delta;
      die $t unless defined $ct->{max};
      die $t unless defined $total_delta;
      if ($total_delta > $ct->{max}) {
        for (my $i = $ct->{max} + 1; $i <= $total_delta; $i++) {
           $ct->{frequencies}->[$i] = 0; 
        }
        $ct->{max} = $total_delta;
      }
      $ct->{frequencies}->[$total_delta]++;
      $ct->{sqtotal} += $total_delta * $total_delta;
    }
  } 
}

sub calculate_derivables {
  my ($self,$group) = @_;
  foreach my $v (values %$group) {
    $self->calculate_cumulative($v);
    $self->calculate_percentiles($v);
    $self->calculate_mean_sd($v);
  }
}

sub calculate_mean_sd {
  my ($self,$counter) = @_;
  my $count = $counter->{count};
  my $total = $counter->{total};
  my $sqtotal = $counter->{sqtotal};
  my $mean = ($count == 0) ? "NaN" : $total / $count;
  my $var = ($count == 0) ? "NaN" : abs($sqtotal/$count - $mean * $mean);
  $var = $var * $count / ($count - 1) if ($count > 0);
  my $sd = sqrt($var);
  $counter->{mean} = sprintf("%1.4f",$mean);
  $counter->{variance} = sprintf("%1.4f",$var);
  $counter->{sd} = sprintf("%1.4f",$sd);
}

sub calculate_cumulative {
  my ($self,$counter) = @_;
  my @frequencies = @{$counter->{frequencies}};
  my $count = $counter->{count};
  if ($count == 0) {
    $counter->{cumulative} = [1];
    return;
  }
  my @result = ();
  my $accu = 0;
  for (my $i = 0; $i < scalar(@frequencies); $i++) {
    $accu += $frequencies[$i];
    $result[$i] = sprintf("%1.4f",$accu / $count);
  }
  $counter->{cumulative} = \@result;
}

sub calculate_percentiles {
  my ($self,$counter) = @_;
  my $count = $counter->{count};
  my @indexes = map { ( $count * $_) / 100 } @{$self->percentiles};
  my @frequencies = @{$counter->{frequencies}};
  my $next = 0; my $last = 0; my $left_count = 0;
  my @values = ();
  while ($next < scalar(@frequencies) && $frequencies[$next] == 0) { $last = $next++ };
  for (my $i = 0; $i < scalar(@indexes); $i++) {
     my $index = $indexes[$i];
     my $ceil_index = ceil($index); 
     my $floor_index = floor($index);
     my $floor_v = $self->skip_percentile(\$left_count,\$floor_index,\$last,\$next,\@frequencies);
     my $ceil_v = $self->skip_percentile(\$left_count,\$ceil_index,\$last,\$next,\@frequencies);
     my $value = $floor_v * (1 - $index + $floor_index) + $ceil_v * (1 - $ceil_index + $index);
     push @values, sprintf("%1.4f",$value);
     $counter->{mean} = $value if ($index >= 0.5 && !exists($counter->{mean})); 
  }
  $counter->{percentiles} = {percents => $self->percentiles,
                            values => \@values};
}

sub skip_percentile {
  my ($self,$left_count,$target,$last,$next,$frequencies) = @_;
  while ($$left_count < $$target && $$next < scalar(@$frequencies)) { 
    if ($frequencies->[$$next] == 0) {
      $$next++;
    }
    else {
      $$last = $$next;
      $frequencies->[$$next]--;
      $$left_count++;
    }
  }
  return $$last;
}

sub gff_iterator {
  my $self = shift;
  my $gff = $self->gff || '/dev/null';
  my $rt = ref($gff);
  return MalariaGEN::AGV::GffUtils::Iterator->new(file => $gff) unless $rt;
  if (blessed($gff)) {
    return $gff if $gff->isa('MalariaGEN::AGV::GffUtils::Iterator');
  }
  return MalariaGEN::AGV::GffUtils::Iterator->new(fh => $gff);
}

sub is_left_out {
  my ($self,$pf) = @_;
  return 0 unless defined $pf && $#$pf >= 0 && ($#$pf > 0 || ($$pf[0] ne "PASS" && $$pf[0] ne "." && $$pf[0] ne "pass"));
  return 1 if $self->filter_all;
  my $mf = $self->filters;
  return 0 unless scalar(@$mf) > 0;
  my %pf = map { $_ => 1 } @$pf;
  foreach my $f (@$mf) {
    exists($pf{$f}) and return 1;
  } 
  return 0;
}

##
#  $object->analyze(alignments => { name => file, name => file})
#
sub analyze {
  my ($self,%args) = @_;
  my $ref = $self->reference;
  my $ref_file = $ref->sequence_file;
  my $alignments = $args{alignments} || $args{alignment} or die "you need to specify some alignment to analyze";
  my $input_pipe = IO::Pipe->new();
  my @names = keys %$alignments;
  my @files = map { $$alignments{$_} } @names;
  $input_pipe->reader("samtools mpileup -f $ref_file " . join(" ",values %$alignments));
  return $self->_analyze_samtools_mpileup(fh => $input_pipe,sample_names => \@names);
}

sub _analyze_samtools_mpileup {
  my ($self,%options) = @_;
  my $input_fh = $options{fh};
  my @sample_names = @{$options{sample_names}};
  my $sample_count = $#sample_names + 1;
  my %all_samples = %{$self->zero_sample_counter_group()};
  my %per_sample = $options{sample_counters} ? map { $_ => $self->zero_sample_counter_group() } @sample_names : ();
  my $gff_iterator = $self->gff_iterator();
  while (my $line = <$input_fh>) {
   chomp $line;
   my ($chr,$pos,$ref,@sample_data) = split(/\t/,$line);
   next unless $ref =~ /[atgc]/i; # N X and so does not count
   my @counters = ('all',$chr);
   push @counters , 'nuclear' unless $chr eq "MT";
   my @types = qw(all);
   my $found = $gff_iterator->seek_pos($chr,$pos);
   push @types, 'featureless' unless $found;
   my $features = $gff_iterator->current_features();
   my $in_gene = 0; my $in_exon = 0;
   foreach my $f (@$features) {
     my $type = $f->{type};
     $in_gene = 1 if ($type eq "gene");
     $in_exon = 1 if ($type eq "exon");
     $in_exon = 1 if ($type eq "CDS");
     $in_exon = 1 if ($type eq "cds");
   }
   push @types, $in_gene ? 'genic' : 'intergenic';
   push @types, $in_exon ? 'exonic' : 'non_coding';
   push @types, 'genic_non_coding' if $in_gene && !$in_exon;
   my $all_coverage = 0;
   for (my $i = 0; $i < $sample_count; $i++) {
     my $sample_name = $sample_names[$i];
     my ($coverage,$pileup,$qual) = @sample_data[($i * 4) .. (($i+1) * 4 - 1)];
     $self->increase_counters(1,$coverage,$per_sample{$sample_name},\@types,\@counters) if exists($per_sample{$sample_name});
     $all_coverage += $coverage;
   }
   $self->increase_counters(1,$all_coverage,\%all_samples,\@types,\@counters);
  }
  foreach my $s (keys %per_sample) {
    foreach my $k (keys %{$per_sample{$s}}) {
      $self->calculate_derivables($per_sample{$s}->{$k});
    }
  }
  foreach my $k (keys %all_samples) {
    $self->calculate_derivables($all_samples{$k});
  }
  return wantarray ? ( all => \%all_samples, %per_sample ) : { all => \%all_samples, %per_sample };
}

#################
# JSON encoding
##

use IO::File;
our $json_coder = undef;

sub write_counter_group {
  $json_coder = $json_coder ||  JSON::XS->new->ascii->pretty->allow_nonref;
  my ($class,$group,$dev) = @_;
  $dev = $dev || '/dev/null';
  my $fh = ref($dev) ? $dev : IO::File->new($dev,"w") or die "could not create file '$dev' to output counter group";
  print $fh $json_coder->encode($group);
  $fh->close();
}

sub read_counter_group {
  my $coder =  JSON::XS->new->ascii->pretty->allow_nonref;
  my ($class,$group,$dev) = @_;
  $dev = $dev || '/dev/null';
 
  my $fh = ref($dev) ? $dev : IO::File->new($dev,"r");
  my $result;
  while (my $line = <$fh>) {
    $coder->incr_parse($line);
  }
  $result = $coder->incr_parse();
  $fh->close();
  return $result;
}

#################
# GD library based plotting
## TODO: must move to some other module:

our $width = 1280 ;
our $height = 1024 ;
our $top = 10 ;
our $bottom = 20 ;
our $left = 10 ;
our $right = 50 ;
our $w = $width - $left - $right ;
our $h = $height - $top - $bottom ;
  
our $maxy = 0 ;
our $max = undef;

sub coverage_plot {  
  my ($self,$counter,$file_name) = @_;
  eval {  
  $max = $counter->{all}->{max} || 1;
  my @c_all = @{$counter->{all}->{frequencies}};
  my @c_coding = @{$counter->{exonic}->{frequencies}};
  my @c_non = @{$counter->{non_coding}->{frequencies}};
  require GD;
  my $im = new GD::Image($width,$height);
  my $white = $im->colorAllocate(255,255,255);
  my $black = $im->colorAllocate(0,0,0);       
  my $red = $im->colorAllocate(200,0,0); 
  my $green = $im->colorAllocate(0,200,0); 
  my $blue = $im->colorAllocate(0,0,200); 

foreach ( 0 .. $max ) {
	$maxy = $c_all[$_] if ( $c_all[$_] || 0 ) > $maxy ;
}
$maxy = 15000 ;

foreach my $cov ( 1 .. $max ) {
	plot_point ( $im,  $cov , $c_all[$cov] || 0 , $red ) ;
	plot_point ( $im, $cov , $c_coding[$cov] || 0 , $green ) ;
	plot_point ( $im, $cov , $c_non[$cov] || 0 , $blue ) ;
}
my $s ;
for ( $s = 1 ; $s < $max ; $s += 1000 ) {
	my $x = getx ( $s ) ;
	$im->line ( $x , $height - $bottom , $x , $height - $bottom + 5 , $black ) ;
	$im->string ( gdSmallFont() , $x , $height - $bottom + 5 , $s , $black ) ;
}
for ( $s = 0 ; $s <= $maxy ; $s += 2500 ) {
	my $x = $width - $right ;
	my $y = $height - $bottom - $s * $h / $maxy ;
	$im->line ( $x , $y , $x + 5 , $y , $black ) ;
	$im->string ( gdSmallFont() , $x , $y , $s , $black ) ;
}

open IMAGEFILE , ">$file_name" ;
binmode IMAGEFILE;
print IMAGEFILE $im->png () ;
close IMAGEFILE ;
};
print STDERR "plot error $file_name: $@\n" if $@;
}

sub iqr {
        my ( $ar , $total , $perc ) = @_ ;
        my $sum = 0 ;
        my $ret = 0 ;
        while ( $sum < $total * $perc / 100 ) {
                $sum += ( $ar->[$ret] || 0 ) ;
                $ret++ ;
        }
        return $ret - 1 ;
}

sub plot_point {
        my ( $im, $cov , $cnt , $col ) = @_ ;
#       return if $cnt <= 0 ;
#       my $l = log ( $cnt ) / log ( 10 ) ;
        my $x = getx ( $cov ) ;
        my $y = $height - $bottom - $cnt * $h / $maxy ;
        $im->setPixel ( $x , $y , $col ) ;
}

sub getx {
        my ( $cov ) = @_ ;
        return $left + log ( $cov ) * $w / log ( $max ) ;
}

1;
