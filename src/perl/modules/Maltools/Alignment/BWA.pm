package Maltools::Alignment::BWA;

use Scalar::Util qw(blessed);
use Moose::Exporter;

use Maltools::Alignment;
use Maltools::Reference;
use File::Temp qw(tempfile);


Moose::Exporter->setup_import_methods(
  as_is => [ 'bwa' ],
);

sub bwa {
  my %args = @_;
  my $in = $args{in} or die "need to indicate an input alignment (or file)";
  my $out = $args{out} or die "need to indicate an output alignment (or file)";
  my $ref = $args{ref} or die "need to indicate a reference";
  my $cpus = $args{cpus} || 1;
  my $sort = $args{sort} || 0;
  my $index = $args{index} || 0;

  $ref = Maltools::Reference->new(sequences => $ref) unless blessed($ref);
  $in = Maltools::Alignment->new(file => $in, reference => $ref) unless blessed($in);
  $out = Maltools::Alignment->new(file => $out, reference => $ref) unless blessed($out); 

  my $rgs = $in->read_groups;
  if (scalar(@$rgs) == 0) {
    
  }

  my @rgs_out = ();


  my $reference = $ref->sequence_file;

  foreach my $rg (@$rgs) {
    my $forward = $rg->forward_fastq_file;
    my $reverse = $rg->reverse_fastq_file;
    my $forward_sai = bwa_aln($forward,$reference,$cpus);
    my $reverse_sai = bwa_aln($reverse,$reference,$cpus);
    my $bwa_sampe =  "bwa sampe $reference $forward_sai $reverse_sai $forward $reverse";
    my $rg_out = $forward;
    my $header_file = create_header_file($ref,$sort,$rg);
    $rg_out =~ s/_1\.fastq$/.bam/;
    my $sorting = $sort ? "picard SortSam INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 SORT_ORDER=coordinate" : "cat";
    my $regrouper_options = regrouper_options($in);
    my $regrouper = "picard AddOrReplaceReadGroups COMPRESSION_LEVEL=0 INPUT=/dev/stdin OUTPUT=/dev/stdout $regrouper_options";
    my $reheader = "picard ReplaceSamHeader INPUT=/dev/stdin HEADER=$header_file OUTPUT=/dev/stdout";
    print STDERR "COMMAND: $bwa_sampe | samtools view -bS - | $sorting | $reheader | $regrouper > $rg_out";

    `$bwa_sampe | samtools view -bS - | $sorting | $reheader | $regrouper > $rg_out`;
    die "could not realign group " . $rg->id  if $?; 
    unlink $forward_sai;
    unlink $reverse_sai; 
    push @rgs_out, $rg_out;
  }  
  ## merge groups 
  my $output = $out->file;
  if (scalar(@rgs_out) == 1) {
    print STDERR "cp $rgs_out[0] $output\n"; 
    `cp $rgs_out[0] $output`; 
    $? and die "failed moving file $rgs_out[0] to $output (code=$?)";
    if ($index) {
      my $output_bai = $output . ".bai";
      `picard BuildBamIndex INPUT=$output OUTPUT=$output_bai`;
      die "could not create index file '$output_bai'" if $?;
    }
  }
  else {
    my $merge_files = join(" ", map { "INPUT=$_" } @rgs_out );
    my $assume_sorted = $sort ? 'true' : 'false';
    my $use_threading = $cpus > 1 ? 'true' : 'false';
    my $sort_order = $sort ? 'coordinate' : 'unsorted';
    my $index_str = $index ? 'true' : 'false';
    `picard MergeSamFiles MERGE_SEQUENCE_DICTIONARIES=true CREATE_INDEX=$index_str USE_THREADING=$use_threading ASSUME_SORTED=$assume_sorted SORT_ORDER=$sort_order $merge_files OUTPUT=$output`;
    die "could not merge files (code=$?)" if $?;
    if ($index) {
      my $output_bai = $output; $output_bai =~ s/.bam$/.bai/;
      `mv $output_bai $output.bai` if -e $output_bai;
    }
  }
}

sub create_header_file {
  my $ref = shift;
  my $sort = shift || 0;
  my $rg = shift;
  my $dict = $ref->dictionary;
  $dict = join("\n", map { "\@SQ\tSN:$$_{name}\tLN:$$_{length}\tM5:$$_{md5}" } @$dict );
  my ($fh,$fn) = tempfile("BwaHdrXXXX", UNLINK => 1);
  my $sort_str = $sort ? "coordinate" : "unsorted";
  my $ref_file = $ref->sequence_file;
  my $rg_id = $rg->id;
  my $rg_line = $rg->line;
  unless ($rg_line =~ s/\tPG:([^\t]+)/PG:bwa_sampe_${rg_id}/) {
    $rg_line .= "\tPG:bwa_sampe_${rg_id}";
  }
  my $library = $rg->library;
  my $platform = $rg->platform;
  my $platform_unit = $rg->platform_unit;
  my $mess = <<EOM
\@HD	VN:1.1	SO:$sort_str
$dict
\@PG	ID:bwa_aln_${rg_id}_1	PN:bwa	CL:bwa aln reference.fa forward.fastq  > forward.sai
\@PG	ID:bwa_aln_${rg_id}_2	PN:bwa	CL:bwa aln reference.fa reverse.fastq  > reverse.sai	PP:bwa_aln_${rg_id}_1
\@PG	ID:bwa_sampe_${rg_id}	PN:bwa	CL:bwa sampe reference.fa forward.sai reverse.sai forward.fastq reverse.fastq | samtools view -bS -	PP:bwa_aln_${rg_id}_2
$rg_line
EOM
 ;
 print $fh $mess,"\n";
 $fh->close(); 
 return $fn;
}

sub bwa_aln {
  my $fastq = shift;
  my $ref = shift;
  my $cpus = shift;
  my $sai = $fastq;
  $sai =~ s/\.fastq$/\.sai/;
  my $cpu_option = $cpus > 1 ? "-t $cpus" : "";
  `bwa aln $cpu_option $ref $fastq > $sai`;
  die "error aligning '$fastq' using reference '$ref'" if $?;
  return $sai;
}

sub underscore {
   my $str = shift;
   return undef unless defined $str;
   $str =~ s/[\(\) \t]+/_/g;
   return $str;
}

sub regrouper_options {
   my $in = shift;
   my @rgs = @{$in->read_groups};
   my $rg = $#rgs >= 0 ? $rgs[0] : Maltools::Alignment::ReadGroups->new();
   my $rg_id = underscore($rg->id || 'UNKOWN');
   my $rg_lb = underscore($rg->library || 'UNKNOWN');
   my $rg_pl = underscore($rg->platform || 'UNKNOWN');
   my $rg_pu = underscore($rg->platform_unit || 'UNKNOWN');
   my $rg_sm = underscore($rg->sample || 'UNKNOWN'); 
   my $rg_cn = underscore($rg->sequencing_center || 'UNKNOWN'); 
   my $rg_ds = underscore($rg->description || 'UNKNOWN');
   my @options = ();
   push @options, "'RGID=$rg_id'" if $rg_id;
   push @options, "'RGPL=$rg_pl'" if $rg_pl;
   push @options, "'RGPU=$rg_pu'" if $rg_pu;
   push @options, "'RGSM=$rg_sm'" if $rg_sm;
   push @options, "'RGCN=$rg_cn'" if $rg_cn;
   push @options, "'RGDS=$rg_ds'" if $rg_ds;
   push @options, "'RGLB=$rg_lb'" if $rg_lb;
  
   my $options = join(" ",@options);
   return $options;
}

1;

