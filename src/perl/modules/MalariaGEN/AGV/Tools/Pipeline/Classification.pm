package MalariaGEN::AGV::Tools::Pipeline::Classification;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift) }
  
our $INPUTS = {
    in_dir => {  type => 'directory', mandatory => 1, mode => 'inout' }, # where to find the genotyping information as generated using the Genotyping stage tool.
# Commented out. For now we relay on the default values set at each tool independently
# TODO add the options to the corresponding commands and make the connections in this tool script code under --DATA--
#    qual_thresholds => { type => 'num', multiple => 1, default => [10,30,100,300] },
#    vc_thr => { type => 'num', default => 30 }, # variant credible threshold, 0 or less means to use the genotyper 'lowqual' filter.
#    vt_qthr  => { type => 'num', default => 20 }, # variant typable threshold, minimum quality applicable.i
#    vt_qthr_hetero = > { type => 'num', default => undef }, # vt_qthr specific for heterogenous calls.
#    vt_qthr_homo => { type => ' num', default => undef }, # vt_qthr specific for homozygous calls.
#    vt_sthr => { type => 'num', default => 0.5 }, # franction of good samples for which the variant must be typable (GQ >= vt_qthr) for consider a variant typable.
#    s_vc_thr => { type => 'num', default => 0.8 }, # fraction of credible variants that the sample should have quality typing (GQ >= vt_qthr) in order to consider the sample good.
};

our $OUTPUTS = { 
    out_dir => { type => 'directory', mandatory => 1 } # where to generat the classification information.
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

# Algorithm outline:
#  1. obtain list of credible variants. (filter by: vc_thr). => cv_list
#  2. obtain list of good samples (thus also bad samples) (using cv_list filter by vt_thr) => gs_list
#  3. obtain list of typable variants (using gs_list and vt_thr filter by vt_sthr) => tv_list
#  4. create new VCF removing poor samples and adding variant classficiation info. 
#  5. do variant bootstraping and distance matrix generataion and make sample trees.
#  6. Repeat 5 just with MT variants (probably better to reconstruct phylogenies)
#  
#  Done along stages 1 to 3:
#  a. generate data for create classification graphs in R. 
#
# Algorithm details:
# 
# 1. obtain list of credible variants (filter by vc_thr) => cv_list.
#  this is virtually in the vcf genotypes file and unless the number of credible variants is some factors less
#  the total number of variants there is no much point in generate a separate list file as VCF offers the functionality.
#  so this step is virtually skipped.
# 
# 2. obtain list of good samples (thus also bad samples)
#  for this it actually make sense to generate a separate list
#  good-sample.list (temporal file)
#  Basically it simply will use the quality counters generated during genotyping and reading the vcf is not 
#  need for that. 
#  
#  Command: agv genotype-query good-samples  -i genotypes.vcf -o good-sample.list -c vc_qthr -t vt_qthr -s svc_thr
# 
# 3 and 4 can be done on the fly together.
# 
#  Command: agv classification -i genotypes.vcf -S good-sample.list -c vc_thr -q vt_qthr -r vt_sthr -o classification.vcf
#
#  we reserve -s for the possibility of indicating individual samples in the command line.
#
# 5. This already programmed using a set of commands:
#
#  agv variant-bootstrap -i classification.vcf -r 100 -s 10000 -o bootstrap_dir 
#  agv variant-distance ...
#  agv sample-tree ...
#
# a. In fact this is quite independent from steps 1 to 3 and should be done directly from 
# the original vcf.  The samples oriented plot can be obtained from the quality counters whereas the
# variant base should be generated from the vcf.

#
# agv sample-cuality-plot -i coverage

# For that we can extend 'classification' command to generate 

# agv classification -i genotypes.vcf -S good-sample.list -c vc_thr -q vt_qthr -r vt_sthr -o classification.vcf -V 10 -V 30 -V 100 -V 300 
# where -V indicate that it should generate information per sample that gives the number of samples at the different vt_thresholds
# provided in 


.PHONY: all

GENOTYPE_DIR={$J->input("in_dir")}
OUTPUT_DIR={$J->output("out_dir") || $J->input("in_dir")}
GENOTYPE_FILE=$(GENOTYPE_DIR)/genotypes.vcf
CLASSIFIED_FILE=$(OUTPUT_DIR)/classified.vcf
CREDIBLE_FILE=$(OUTPUT_DIR)/credible.vcf
TYPABLE_FILE=$(OUTPUT_DIR)/typable.vcf
SAMPLE_TABLE_FILE=$(OUTPUT_DIR)/sample-table.tab
GOOD_SAMPLES_FILE=$(OUTPUT_DIR)/good-samples.list

# Global
BOOTSTRAP_DIR=$(OUTPUT_DIR)/bootstrap
BOOTSTRAP_SUMMARY_FILE=$(BOOTSTRAP_DIR)/bootstrap-summary.txt
BOOTSTRAP_FILE_LIST=$(BOOTSTRAP_DIR)/bootstrap-file.list
DISTANCE_FILE=$(OUTPUT_DIR)/sample-distance.txt
CONSENSUS_TREE_FILE=$(OUTPUT_DIR)/consensus.newick
NJ_TREES_FILE=$(OUTPUT_DIR)/nj-trees.newick
TREE_FILE=$(OUTPUT_DIR)/sample-tree.xml

# MT only 
MT_BOOTSTRAP_DIR=$(OUTPUT_DIR)/bootstrap_mt
MT_BOOTSTRAP_SUMMARY_FILE=$(MT_BOOTSTRAP_DIR)/bootstrap-summary.txt
MT_BOOTSTRAP_FILE_LIST=$(MT_BOOTSTRAP_DIR)/bootstrap-file.list
MT_DISTANCE_FILE=$(OUTPUT_DIR)/mt-sample-distance.txt
MT_CONSENSUS_TREE_FILE=$(OUTPUT_DIR)/mt-consensus.newick
MT_NJ_TREES_FILE=$(OUTPUT_DIR)/mt-nj-trees.newick
MT_TREE_FILE=$(OUTPUT_DIR)/mt-sample-tree.xml

all: $(SAMPLE_TABLE_FILE) $(TYPABLE_FILE) $(TREE_FILE) $(MT_TREE_FILE)

$(SAMPLE_TABLE_FILE): $(GENOTYPE_FILE)
	agv sample-table -i $< -o $@

$(GOOD_SAMPLES_FILE): $(SAMPLE_TABLE_FILE)
	cat $< | grep -v "^#" | awk '($$2 = "GOOD")\{ print $$1 \}' > $@

$(CLASSIFIED_FILE): $(GOOD_SAMPLES_FILE) $(GENOTYPE_FILE)
	agv variant-info -i $(GENOTYPE_FILE) CLASS:gs_file=$< SCQ:gs_file=$<:thr=10,20,30,40,50,60 SCD:gs_file=$<:10,20,30,40,50,60 -o $@

$(CREDIBLE_FILE): $(CLASSIFIED_FILE)
	agv variant-slice -i $< "/defined(CLASS) && (CLASS eq 'Credible' || CLASS eq 'Typable')/+" -o $@

$(TYPABLE_FILE): $(CREDIBLE_FILE) $(GOOD_SAMPLES_FILE)
	agv variant-slice -i $< -S $(GOOD_SAMPLES_FILE) "/CLASS eq 'Typable'/+" -o $@

$(BOOTSTRAP_SUMMARY_FILE): $(TYPABLE_FILE)
	mkdir -p $(BOOTSTRAP_DIR)
	agv-local variant-bootstrap -i $< -o $(BOOTSTRAP_DIR) -S $@ -r 100 -s 1000 --seed 13 

$(BOOTSTRAP_FILE_LIST) : $(BOOTSTRAP_SUMMARY_FILE)
	find $(BOOTSTRAP_DIR) -name "*.vcf" > $@

$(DISTANCE_FILE) : $(BOOTSTRAP_FILE_LIST)
	agv-local sample-distance -I $< -o $@

$(CONSENSUS_TREE_FILE) : $(DISTANCE_FILE) 
	agv-local sample-tree -i $< -C $@ -T $(NJ_TREES_FILE) -c 

$(TREE_FILE) : $(CONSENSUS_TREE_FILE) 
	agv-local annotate-sample-tree -i $< -o $@ -I newick -O phyloxml -d property

$(MT_BOOTSTRAP_SUMMARY_FILE): $(TYPABLE_FILE)
	mkdir -p $(MT_BOOTSTRAP_DIR)
	agv-local variant-bootstrap -i $< -o $(MT_BOOTSTRAP_DIR) -S $@ -c MT -r 100 -s 1.0 --seed 13

$(MT_BOOTSTRAP_FILE_LIST) : $(MT_BOOTSTRAP_SUMMARY_FILE)
	find $(MT_BOOTSTRAP_DIR) -name "*.vcf" > $@

$(MT_DISTANCE_FILE) : $(MT_BOOTSTRAP_FILE_LIST)
	agv-local sample-distance -I $< -o $@

$(MT_CONSENSUS_TREE_FILE) : $(MT_DISTANCE_FILE)
	agv-local sample-tree -i $< -C $@ -T $(MT_NJ_TREES_FILE) -c

$(MT_TREE_FILE) : $(MT_CONSENSUS_TREE_FILE)
	agv-local annotate-sample-tree -i $< -o $@ -I newick -O phyloxml -d property

