package MalariaGEN::AGV::Tools::Gatk::VariantFiltration;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;
use POSIX;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) }
  
our $INPUTS = {
    #TODO revert obvious problematic changes
    in => { type => 'vcf', multiple => 1, mandatory => 1 }, 
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
  $args{cpu_time} = $self->calulcate_cpu_time(%args) unless exists $args{cpu_time};
  return $self->SUPER::job(%args);
}

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  my $sample_count = scalar(@{$args{inputs} || []});
  my $cpus = 1 + floor(log($sample_count) / log(1.5));
  $cpus = $cpus + 1 unless $cpus % 2 == 0;
  $cpus = 12 if $cpus > 12;
  return $cpus;
}

sub calculate_memory {
  my ($self,%args) = @_;
  my $cpu_ratio = $args{cpu_ratio} || $self->calculate_cpu_ratio(%args);
  return ceil(1000 * $cpu_ratio);
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  my $sample_count = scalar(@{$args{inputs} || [1]});
  # very adhoc need to actually give a better empirical estimate for these figures also based in reference sequence length and BAM file size.
  my $days = floor(1 + $sample_count / 10); 
  return 60 * 60 * 24 * $days;
}

sub script_template {
   my ($self,%args) = @_;
   return TTS('gatk --memory {$J->memory} -T UnifiedGenotyper -R {$J->input("ref")} { join(" ",map {"-I $_"} @{$J->input("samples")} ); } -o {$J->output("vcf")}');
}

1;
