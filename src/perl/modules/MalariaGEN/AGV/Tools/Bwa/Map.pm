package MalariaGEN::AGV::Tools::Bwa::Map;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use Text::Template;

our $INPUTS = {
    in => { type => 'bam', indexed => 0 },
    ref => {  type => { name => 'fasta', indexed => 1 }, mandatory => 1 },
};

our $OUTPUTS = { 
    out => { type => { 'bam', indexed => 1 }, mandatory => 1 }
};

sub interpreter {
   return '$MAKE';
}

1;

__DATA__
'{
  $in_bam = $J->input("in");
  $f_sai = $E->tempfile();
  $r_sai = $E->tempfile();
  $tmp1 = $E->tempfile("aln_XXXXX",".sam");
  $f_fq = $J->tempfile("f_fq_XXXXX",".fastq");  
  $r_fq = $J->tempfile("r_fq_XXXXX",".fastq"); 
  $fq_ts = $J->tempfile();
  $db = $J->input("ref");
  $db_fai = $J->input("ref") ".fai";
  unless (defined $db_fai) { $db_fai = $db; $db_fai =~ s/(\.((fa)|(fasta)))?(\.gz)?$/.fai/; } 
  $out= $J->output("out"); 
  $bam_out = $out =~ /\.bam$/ ? $out : $J->tempfile("bam_XXXXX",".bam"); "";
}.PHONY : rule_0 rule_1 rule_2 

{ $f_fq } { $r_fq } : { $fq_ts }

{ $fq_ts }: { $in_bam }
	agv run picard SamToFastq INPUT=$< FASTQ={$f_fq} SECOND_END_FASTQ={$r_fq} 
	date > $@

rule_0: rule_1 rule_2 {$db_fai}
	agv run bwa sampe {$db} {$f_sai} {$r_sai} {$f_fq} {$r_fq} > {$tmp1}
	agv run picard SortSam INPUT={$tmp1} OUTPUT={$bam_out} CREATE_INDEX=true
        { $out eq $bam_out  ? "" : "mv $bam_out $out; mv ${bam_out}.bai ${out}.bai" }
	rm $tmp1	

rule_1: {$f_fq}
	bwa aln {$db} {$f_fq} > {$f_sai}

rule_2: {$r_fq}
	bwa aln {$db} {$r_fq} > {$r_sai}

