package MalariaGEN::AGV::Tools::Gatk::Genotyper;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;
use POSIX;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) }
  
our $INPUTS = {
    genotyper => { type => 'string', default => 'UnifiedGenotyper' },
    reference => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
    samples => { type => { name => 'bam', indexed => 1}, multiple => 1, mandatory => 1 },
    intervals => { type => 'string', multiple => 1, mandatory => 0, default => []},
    annotations => { type => 'string', multiple => 1, mandatory => 0, default => [],},
    options => { type => 'hash', mandatory => 0, default => {},},
    rods => { type => 'hash', mandatory => 0, default => {},},
};

our $OUTPUTS = { 
    out => { type => 'vcf', mandatory => 1 }
};


sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub job {
  my ($self,%args) = @_;
  $args{cpu_ratio} = $self->calculate_cpu_ratio(%args) unless exists $args{cpu_ratio};
  $args{memory} = $self->calculate_memory(%args) unless exists $args{memory};
  $args{cpu_time} = $self->calculate_cpu_time(%args) unless exists $args{cpu_time};
  return $self->SUPER::job(%args);
}

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1 if (scalar(@{$args{inputs}->{regions} || []}) > 0);
  my $sample_count = scalar(@{$args{inputs}->{samples} || []});
  my $cpus = 1 + floor(log($sample_count) / log(1.5));
  $cpus = $cpus + 1 unless $cpus % 2 == 0;
  $cpus = 12 if $cpus > 12;
  return $cpus;
}

sub calculate_memory {
  my ($self,%args) = @_;
  my $cpu_ratio = $args{cpu_ratio} || $self->calculate_cpu_ratio(%args);
  return ceil(2000 * $cpu_ratio);
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  return 3 * 60 * 60  if (scalar(@{$args{inputs}->{regions} || []}) > 0);
  my $sample_count = scalar(@{$args{inputs}->{samples} || [1]});
  # very adhoc need to actually give a better empirical estimate for these figures also based in reference sequence length and BAM file size.
  my $days = floor(1 + log($sample_count / 10)); 
  return 60 * 60 * 24 * $days;
}

sub script_template {
   my ($self,%args) = @_;
   return TTS('gatk --memory {$J->memory} {$J->cpu_count == 1 ? "":"-nt " . $J->cpu_count . " " }'
     . '-T {$J->input("genotyper")} ' 
     . '-R {$J->input("reference")} ' 
     . ' { join(" ",map {"-A $_"} @{$J->input("annotations")} ); } ' 
     . ' { join(" ",map {"-I $_"} @{$J->input("samples")} ); } '
     . ' { join(" ",map {"-L $_"} @{$J->input("intervals")} ); } ' 
     . ' { join(" ",map {"-$_ " . $J->input("options")->{$_} } keys %{$J->input("options")}); } '
     . '{ $ref_base = $J->input("reference"); $ref_base =~ s/\.fa(sta)?$//; "" }'
     . ' { join(" ",map {"-B:$_," . $J->input("rods")->{$_}->{type} . " " . $ref_base . $J->input("rods")->{$_}->{extension} } keys %{$J->input("rods")} )}'
     . ' -o {$J->output("out")}' . "\n" );
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;
