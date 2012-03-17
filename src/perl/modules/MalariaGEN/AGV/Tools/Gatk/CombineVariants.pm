package MalariaGEN::AGV::Tools::Gatk::CombineVariants;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;

  
our $INPUTS = {
    in => { type => { name => 'vcf' } , multiple => 1, mandatory => 1 },
    ref => { type => { name => 'fasta', indexed => 1 }, mandatory => 1 },
};

our $OUTPUTS = { 
    out => { type => { name => 'vcf', indexed => 1 }, mandatory => 1 }
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
  return 1;
}

sub calculate_memory {
  my ($self,%args) = @_;
  return 1024;
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  return 60 * 60;
}

sub interpreter {
   my ($self) = @_;
   return '/bin/bash';
}

1;

__DATA__
{ 
  $idx = 0;
  @input_arguments = map {"-B:input" . $idx++ . ",VCF $_"} @{$J->input("in")};
  $out = $J->output("out");
  $tmp = $J->tempfile('binds_XXXX');
  open my $tmp_fh , $tmp;
  print $tmp_fh $_,"\n" foreach (@input_arguments);
  close $tmp_fh;
  $ref = $J->input("ref"); ''
}
gatk --memory {$J->memory} -T CombineVariants -R {$ref} {$input_arguments} -o {$out} -setKey null
