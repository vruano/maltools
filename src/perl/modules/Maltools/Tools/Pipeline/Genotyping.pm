package Maltools::Tools::Pipeline::Genotyping;

use base 'Maltools::Tool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift) }
  
our $INPUTS = {
    time => { type => 'string', mandatory => 0 },
    group => {  type => 'string', mandatory => 1 },
    uq_file => { type => 'file', mandatory => 1},
    gff_file => { type => 'gff', mandatory => 1},
    alignment_class => { type => 'string', default => 'realigned' },
};

our $OUTPUTS = { 
    out_dir => { type => 'directory', mandatory => 1 }
};

sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
   return '$MAKE'; 
}

sub script_template {
   return Text::Template->new (TYPE => 'FILEHANDLE', SOURCE => \*DATA);
}

sub command_template {
	   return TTS('{$T->interpreter} -j 2 -f {$S} all');
}

1;

__DATA__

{ $time = $J->input('time') || ''; '' }
.DELETE_ON_ERROR:
.PHONY: all

GROUP={$J->input("group")}
RAW_GENOTYPES={$J->output("out_dir") . "/raw-genotypes.vcf" }
TEMP_FILE={$J->tempfile()}
TEMP_FILE_2={$J->tempfile()}
COVERAGE_DIR={ $J->output("out_dir") . "/coverage" }
UQ_FILE={$J->input("uq_file")}
GFF_FILE={$J->input("gff_file")}
GENOTYPES_FILE={$J->output("out_dir") . "/genotypes.vcf" }
ALIGNMENT_CLASS={$J->input("alignment_class")}

all: $(GENOTYPES_FILE)

SAMPLE_LIST=$(shell agv query group-samples { $time ? "--time $time" : ''} -H $(GROUP)) 

$(GENOTYPE_FILE) : 
	agv genotype --fragments 100 -o $@ -c $(ALIGNMENT_CLASS) $(SAMPLE_LIST)
#	agv genotype $(SAMPLE_LIST) -o $(TEMP_FILE) -c $(ALIGNMENT_CLASS)
#	agv variant-info CODING:gff_file=$(GFF_FILE) UQ:uq_file=$(UQ_FILE) -i $(TEMP_FILE) -o $(TEMP_FILE_2)
#	agv variant-filter BIAL SNPQ mALF:thr=0.01 mALSRF:thr=10 -i $(TEMP_FILE_2) -o $(RAW_GENOTYPES)

#$(COVERAGE_DIR)/cvg-all.json : $(RAW_GENOTYPES)
#	agv coverage-analysis -F -g $(GFF_FILE) -i $(RAW_GENOTYPES) -o $(COVERAGE_DIR)

#$(GENOTYPES_FILE) : $(RAW_GENOTYPES) $(COVERAGE_DIR)/cvg-all.json
#	agv variant-info -i $(RAW_GENOTYPES) CCP:cvg_dir=$(COVERAGE_DIR) |\
#	agv variant-filter -o $@ hiCVG:thr=0.85 loCVG:thr=0.15 
