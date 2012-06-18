package Maltools::Tools::Pipeline::LaneAlign;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;

our $INPUTS = {
    in_irods => { type => 'string', mandatory => 1},
    in => { type => 'bam',  mandatory => 1 },
    ref => {  type => { name => 'reference', indexed => 1, sam_hdr => 1 }, mandatory => 1 },
    pg_id => { type => 'string', mandatory => 1 },
    pg_vn => { type => 'string', mandatory => 1 },
    pg_pn => { type => 'string', mandatory => 1 },
    pg_ds => { type => 'string', mandatory => 1 },
    pg_cl => { type => 'string', mandatory => 1 },
    rg_id => { type => 'string', mandatory => 1},
    rg_lb => { type => 'string', mandatory => 1},
    rg_pl => { type => 'string' },
    rg_pu => { type => 'string' },
    rg_cn => { type => 'string', default => 'Wellcome Trust Sanger Institute' },
    rg_ds => { type => 'string', default => 'Anopheles Genome Variation project produced lane level alignment' },
    rg_sm => { type => 'string', mandatory => 1},
    rg_pu => { type => 'string' },
};

our $OUTPUTS = { 
    out => { type => { name => 'bam', indexed => 1}, mandatory => 1 }
};

sub calculate_cpu_ratio {
   return 4;
}

sub calculate_cpu_time {
   return 4 * 60 * 60 * 4;
}

sub calculate_memory {
   return 2000;
}

sub interpreter {
   return '$SHELL';
}

1;

__DATA__
{ use POSIX;
  $in_bam = $J->input("in");
  $in_irods = $J->input("in_irods");
  $from_irods = $in_irods && ! -e $in_bam;
  $in_bam = $in_bam ? $in_bam : $J->tempfile("bam_XXXXX",".bam");
  $f_sai = $E->tempfile();
  $cpu_count = $J->cpu_count;
  $h_cpu_count = ceil($cpu_count / 2);
  $r_sai = $E->tempfile();
  $tmp1 = $E->tempfile("aln_XXXXX",".sam");
  $f_fq = $J->tempfile("f_fq_XXXXX",".fastq");  
  $r_fq = $J->tempfile("r_fq_XXXXX",".fastq"); 
  $fq_ts = $J->tempfile();
  $db = $J->input("ref");
  $db_fai = $J->input("ref") . ".fai";
  unless (defined $db_fai) { $db_fai = $db; $db_fai =~ s/(\.((fa)|(fasta)))?(\.gz)?$/.fai/; } 
  $out= $J->output("out"); 
  $bam_out = $out =~ /\.bam$/ ? $out : $J->tempfile("bam_XXXXX",".bam");
  $bai_out = $bam_out; $bai_out =~ s/\.bam$/.bai/; 
  $rg_id = $J->input("rg_id") ? "'" . $J->input("rg_id") . "'" : '${RG_ID}';
  $rg_lb = $J->input("rg_lb") ? "'" . $J->input("rg_lb") . "'" : '${RG_LB}';
  $rg_cn = $J->input("rg_cn") ? "'" . $J->input("rg_cn") . "'" : '${RG_CN}';
  $rg_sm = $J->input("rg_sm") ? "'" . $J->input("rg_sm") . "'" : '${RG_SM}';
  $rg_pu = $J->input("rg_pu") ? "'" . $J->input("rg_pu") . "'" : '${RG_PU}';
  $rg_ds = $J->input("rg_ds") ? "'" . $J->input("rg_ds") . "'" : '${RG_DS}';
  $rg_pl = $J->input("rg_pl") ? "'" . $J->input("rg_pl") . "'" : '${RG_PL}';
  $rg_id =~ s/\s/_/g;
  $rg_lb =~ s/\s/_/g;
  $rg_cn =~ s/\s/_/g;
  $rg_sm =~ s/\s/_/g;
  $rg_pu =~ s/\s/_/g;
  $rg_ds =~ s/\s/_/g;
  $rg_pl =~ s/\s/_/g;
  $pg_id = $J->input("pg_id");
  $pg_vn = $J->input("pg_vn");
  $pg_cl = $J->input("pg_cl");
  $pg_ds = $J->input("pg_ds");
  $pg_pn = $J->input("pg_pn");
  $pg_line = join("\t",'@PG',"ID:$pg_id","PN:$pg_pn","VN:$pg_vn","DS:$pg_ds","CL:$pg_cl"); "" }

{$from_irods ? "irsync -N 0 -s 'i:$in_irods' $in_bam" : "" } || exit 1
samtools view -H {$in_bam} > {$tmp1}.old-hdr || exit 2
set RG_PU=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "PU:[^{"\t"}]*" | sed 's/PU://g;s/\s/_/g'`
set RG_PL=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "PL:[^{"\t"}]*" | sed 's/PL://g;s/\s/_/g'`
set RG_ID=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "ID:[^{"\t"}]*" | sed 's/ID://g;s/\s/_/g'`
set RG_LB=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "LB:[^{"\t"}]*" | sed 's/LB://g;s/\s/_/g'`
set RG_CN=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "CN:[^{"\t"}]*" | sed 's/CN://g;s/\s/_/g'`
set RG_SM=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "SM:[^{"\t"}]*" | sed 's/SM://g;s/\s/_/g'`
set RG_DS=`cat {$tmp1}.old-hdr | grep "@RG" | grep -o "DS:[^{"\t"}]*" | sed 's/DS://g;s/\s/_/g'`

picard RevertSam RESTORE_ORIGINAL_QUALITIES=false COMPRESSION_LEVEL=0 INPUT={$in_bam} OUTPUT=/dev/stdout |\
picard SamToFastq INPUT=/dev/stdin FASTQ={$f_fq} SECOND_END_FASTQ={$r_fq} || exit 3;
echo "{$pg_line}" | cat {$db}.sam-hdr - | merge_irods_header {$tmp1}.old-hdr > {$tmp1}.hdr || exit 4;
bwa aln { $cpu_count <= 1 ? '' : "-t $cpu_count" } {$db} {$f_fq}  > {$f_sai} || exit 5;
bwa aln { $cpu_count <= 1 ? '' : "-t $cpu_count" } {$db} {$r_fq}  > {$r_sai} || exit 6;
bwa sampe {$db} {$f_sai} {$r_sai} {$f_fq} {$r_fq} |\
samtools view -bS - |\
picard SortSam COMPRESSION_LEVEL=0 INPUT=/dev/stdin SORT_ORDER=coordinate OUTPUT=/dev/stdout |\
picard ReplaceSamHeader COMPRESSION_LEVEL=0 INPUT=/dev/stdin HEADER={$tmp1}.hdr OUTPUT=/dev/stdout |\
picard AddOrReplaceReadGroups INPUT=/dev/stdin OUTPUT={$bam_out} RGID={$rg_id} RGLB={$rg_lb} RGPL={$rg_pl} RGCN={$rg_cn} RGSM={$rg_sm} RGPU={$rg_pu} RGDS={$rg_ds} CREATE_INDEX=true || exit 7; 

{ $out eq $bam_out  ? "mv $bai_out ${out}.bai || exit 8" : "mv $bam_out $out; mv $bai_out ${out}.bai || exit 8" }
rm { $tmp1 }* {$f_fq} {$r_fq} || exit 9;


