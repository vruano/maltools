package Maltools::Tools::Gatk::CompareAlignments;

use base 'Maltools::Tool';
use strict;
use warnings;
use Text::Template;
use POSIX;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) }
  
our $INPUTS = {
    left => { type => { name => 'bam', indexed => 1 }, mandatory => 1},
    right => { type => { name => 'bam', indexed => 1 }, mandatory => 1 },
    ref => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
};

our $OUTPUTS = { 
    out => { type => 'file', mandatory => 1 }
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
  return 2000;
}

sub calculate_cpu_time {
  my ($self,%args) = @_;
  return 60 * 60;
}

sub script_template {
   my ($self,%args) = @_;
   return TTS('gatk --memory {$J->memory} {$J->cpu_count == 1 ? "":"-nt " . $J->cpu_count . " " }'
     . '-T CompareAlignments ' 
     . '-R {$J->input("ref")} '
     . '-I {$J->input("left")} '
     . '-I {$J->input("right")} '
     . '-o {$J->output("out")}' . "\n" );
}

sub interpreter {
   my ($self) = @_;
   return '$SHELL';
}

1;
