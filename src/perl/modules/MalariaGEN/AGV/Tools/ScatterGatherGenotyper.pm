package MalariaGEN::AGV::Tools::ScatterGatherGenotyper;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;
use POSIX;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) }
  
our $INPUTS = {
    fragments => { type => 'num', default => 100 },
    threads_per_fragment => { type => 'num', default => 2 },
    ref => {  type => { name => 'fasta', indexed => 1 },  mandatory => 1 },
    samples => { type => { name => 'bam', indexed => 1}, multiple => 1, mandatory => 1 },
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
  return 0.1;
}

sub calculate_memory {
  my ($self,%args) = @_;
  return 100;
}

sub calculate_cpu_time {
  return 60 * 60;
}

sub command_template {
  return '{$J->interpreter()} -j {int($J->input("fragments") * 2)} -k {$S}';
} 

sub interpreter {
   return '$MAKE';
}

1;

__DATA__

.DELETE_ON_ERROR:
.PHONY: all

TIME={$J->input("time")}
OUTFILE={$J->output("out")}
FRAGMENTS={$J->input("fragments")}
THREADS={$J->input("threads")}
SAMPLES={join(" ",@{$J->input("samples")})}

ALIGNMENT_CLASS=realigned


INTERVAL_OUT_FILE = $(subst :,_,$(addprefix $(OUTFILE).,$(1)))
INTERVALS := $(shell agv reference intervals -n $(FRAGMENTS))
INTERVAL_OUT_FILES := $(foreach interval,$(INTERVALS),$(call INTERVAL_OUT_FILE,$(interval)))

define INTERVAL_rule

$(call INTERVAL_OUT_FILE,$(1)): 
	agv genotype $$(SAMPLES) -t $$(THREADS) -c $$(ALIGNMENT_CLASS) -r $(1) -o $$@

endef

$(eval $(foreach interval,$(INTERVALS),$(call INTERVAL_rule,$(interval))))

$(info $(foreach interval,$(INTERVALS),$(call INTERVAL_rule,$(interval))))
$(info $(OUTFILE) $(INTERVAL_OUT_FILES))

$(OUTFILE): $(INTERVAL_OUT_FILES)
	vcf-concat.pl $(INTERVAL_OUT_FILES) > $@ 
	rm $^ 

all : $(OUTFILE)
