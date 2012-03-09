package MalariaGEN::AGV::Commands::CrossVcf;

use Moose;
extends 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(sanger_config);
use Getopt::Long qw(:config);

use IO::File;
use DateTime;
use DateTime::Format::ISO8601;
use File::Spec::Functions qw(catfile);
use base 'MalariaGEN::AGV::Command';
use Digest::MD5;
use Scalar::Util qw(blessed);
use Vcf;

sub help_summary {
   return 'transform the results of th cross analysis into a VCF formated file';
}

sub help_text {
   return "cross-vcf:\n\t" .
          "tranform cross error calls into vcf \n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " cross-vcf -m matrix-dir -r report-dir  -o output-vcf\n";
}

sub execute {
  my ($self,$samtrak_site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my ($output_file, $matrix_dir, $report_dir, $original_vcf);
  my @parents = ();
  GetOptions(
     "o=s" => \$output_file,
     "m=s" => \$matrix_dir,
     "r=s" => \$report_dir,
     "v=s" => \$original_vcf,
     "p=s@" => \@parents,
  );

  $#parents == 1 or return $self->error_return("you must specify two parents");

  $matrix_dir or return $self->error_return("you must provide an non-empty matrix directory name");
  $report_dir or return $self->error_return("you must provide an non-empty report directory name");
  $original_vcf or return $self->error_return("you must provide a non-empty original vcf");
  $output_file ||= "-";
  
  -d $matrix_dir or return $self->error_return("input matrix directory does no exists or is not a directory '$matrix_dir'");
  -d $report_dir or return $self->error_return("input report directory does no exists or is not a directory '$report_dir'");
  -f $original_vcf or return $self->error_return("inpu original vcf does no exists or is not a regular file '$original_vcf'"); 

  my $out_fh =  $output_file eq "-" ? \*STDOUT : IO::File->new($output_file,"w") or return $self->error_return("could not open output '$output_file' for writing");
  my $in_vcf = Vcf->new(file => $original_vcf);
  
  $in_vcf->parse_header();
  $in_vcf->add_header_line({key => 'INFO', ID=> 'SEGREGATING', Type => 'Flag' , Number => 0,
      Description => 'Indicates whether the site is segregating or not based on parent genotypes' });
  $in_vcf->add_header_line({key => 'FORMAT', ID => 'ER', Type => 'Integer', Number => 1, 
      Description => 'Indicates whether the call was an error by indicating the type (0=No error, 1=Mendelian error, 2=Errand error)' });
 
  my @samples = $in_vcf->get_samples;
  my $next_index = 0;
  my %sample_to_index = map { $_ => $next_index++ } @samples;
  my %parent_to_index = map { $_ => $sample_to_index{$_} || -1 } @parents;
  my @progeny = grep { ! exists($parent_to_index{$_}) } @samples;
  my %progeny_to_index = map { $_ => $sample_to_index{$_} } @progeny;

  my @missing_parents = grep { $parent_to_index{$_} == -1 } @parents;

  $#missing_parents < 0 or return $self->error_return("there are some missing parent samples not present in the original vcf: "  . join(", ",@missing_parents)); 

  print $out_fh $in_vcf->format_header();
  my $variant = $in_vcf->next_data_hash;
  
  my $cross_pointer = undef;

  while ($variant) {
    $cross_pointer = $self->_cross_data_search($cross_pointer,$matrix_dir,$report_dir,$variant);
    unless ($cross_pointer->{found}) {
       $variant = $in_vcf->next_data_hash;
       next;
    }
    
    my $gtypes = $variant->{gtypes};
    my @parents_gts = split(/,/,$cross_pointer->{parents_line});
    $#parents_gts == $#parents or return $self->error_return("missmatch in the number of parents and parents genotypes at $$variant{CHROM} $$variant{POS}: " 
                 . ($#parents_gts + 1) . " != " . ($#parents + 1));
    my %parent_gts_count = ();
    foreach my $parent (@parents) {
      $gtypes->{$parent}->{GT} = shift @parents_gts;
      $parent_gts_count{$gtypes->{$parent}->{GT}}++;
    }
 
    my $segregating = scalar(keys %parent_gts_count) == 2;
    if ($segregating) { $variant->{INFO}->{SEGREGATING} = undef };
    
    my @progeny_gts = split(/,/,$cross_pointer->{progeny_line});
    my @mendels_err = split(/\t/,$cross_pointer->{mendels_line} || "");
    my @errands_err = split(/\t/,$cross_pointer->{errands_line} || "");
    $#progeny_gts == $#progeny or return $self->error_return("missmatch in the number of progeny and progeny genotypes at $$variant{CHROM} $$variant{POS}");
    foreach my $progeny (@progeny) {
      my $mendel = shift @mendels_err || 0;
      my $errand = shift @errands_err || 0;
      my $gt = shift @progeny_gts;
      my $er = $mendel + 2 * $errand;
      $gt = ($er == 0) ? $gt : ($gt == 0 ? 1 : 0); 
      $gtypes->{$progeny}->{GT} = $gt;
      $gtypes->{$progeny}->{ER} = $er;
    }
    push @{$variant->{FORMAT}},'ER';

    print $out_fh $in_vcf->format_line($variant);
    $variant = $in_vcf->next_data_hash;
  }
  $out_fh->close();
  return $self->ok_return();
 
}

sub _cross_data_search {
  my ($self,$pointer,$matrix_dir,$report_dir,$variant) = @_;

  my $target_chr = $variant->{CHROM};
  unless ($pointer) {
     $pointer = { 
         chr => $target_chr,
         pos => 0,
         positions_fh => IO::File->new(catfile($matrix_dir,"positions.$target_chr"),"r") || IO::File->new("/dev/null"),
         parents_fh => IO::File->new(catfile($matrix_dir,"parents.$target_chr"),"r") || IO::File->new("/dev/null"),
         progeny_fh => IO::File->new(catfile($matrix_dir,"progeny.$target_chr"),"r") || IO::File->new("/dev/null"),
         errands_fh => IO::File->new(catfile($report_dir,"errant.${target_chr}.txt"),"r") || IO::File->new("/dev/null"),
         mendels_fh => IO::File->new(catfile($report_dir,"mendels.${target_chr}.txt"),"r") || IO::File->new("/dev/null"), 
     };
  }
  elsif ($pointer->{chr} ne $target_chr) {
    $pointer->{$_}->close() foreach (qw(parents_fh progeny_fh errands_fh mendels_fh));
    $pointer->{pos} = 0;
    $pointer->{chr} = $target_chr;
    $pointer->{positions_fh} = IO::File->new(catfile($matrix_dir,"positions.$target_chr"),"r") || IO::File->new("/dev/null");
    $pointer->{parents_fh} = IO::File->new(catfile($matrix_dir,"parents.$target_chr"),"r") || IO::File->new("/dev/null");
    $pointer->{progeny_fh} = IO::File->new(catfile($matrix_dir,"progeny.$target_chr"),"r") || IO::File->new("/dev/null");
    $pointer->{errands_fh} = IO::File->new(catfile($report_dir,"errant.${target_chr}.txt"),"r") || IO::File->new("/dev/null");
    $pointer->{mendels_fh} = IO::File->new(catfile($report_dir,"mendels.${target_chr}.txt"),"r") || IO::File->new("/dev/null");
  }

  my ($positions_fh,$parents_fh,$progeny_fh,$errands_fh,$mendels_fh) = 
     map { $pointer->{$_} } (qw(positions_fh parents_fh progeny_fh errands_fh mendels_fh)); 

  while ($pointer->{pos} < $variant->{POS}) {
    $pointer->{positions_line} = <$positions_fh> || 'inf';
    $pointer->{parents_line} = <$parents_fh>;
    $pointer->{progeny_line} = <$progeny_fh>;
    $pointer->{errands_line} = <$errands_fh>;
    $pointer->{mendels_line} = <$mendels_fh>;
    foreach (qw(positions_line progeny_line parents_line errands_line mendels_line)) {
      chomp $pointer->{$_} if defined $pointer->{$_};
    }
    $pointer->{pos} = $pointer->{positions_line};
  }
  $pointer->{found} = $pointer->{pos} == $variant->{POS};
 
  return $pointer;
}

1;
