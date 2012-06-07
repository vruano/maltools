package Maltools::Tools::Gatk::CombineCounts;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;
use POSIX;
  
our $INPUTS = {
    in => { type => { name => 'file' } , multiple => 1, mandatory => 1 },
};

our $OUTPUTS = { 
    out => { type => { name => 'file' }, mandatory => 1 }
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
  return 500;
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  return 5 * 60;
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;

__DATA__
{ 
  $input_arguments = join(" ",map { "-I $_" } @{$J->input("in")});
  $out = $J->output("out"); '' };
gatk --memory {$J->memory} --main net.malariagen.gatk.math.CombineIntegerDistributionSets {$input_arguments} -o {$out}
