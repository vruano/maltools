package Maltools::Commands::Genotype;

use base 'Maltools::Command';

use strict;
use warnings;

use Maltools::Config qw(data_config 
       sanger_config jobs_config);
use Maltools::Tool;
use Getopt::Long;
use IO::File;

my %genotypers = (
      ReadCounts => {tool => 'gatk-sgReadCounter'},
      Samtools => 1,
      SamtoolsLegacy => 1,
      DiploidGenotyper => {tool => 'gatk-sgGenotyper'},
      MetaGenotyper => {
         tool => 'gatk-genotyper',
         genotyper => 'MetaGenotyper',
         annotations => [qw(
           AverageBaseQuality
           ReadDepthAndAllelicFractionBySample
           DepthPerAlleleByVariant
           UniquenessScore
           NumberSamplesWithData)],
         options => {
           mmq => 20,
           mbq => 20,
           smodel => 'Bernoulli/0.0001/',
           mgc => -200000,
           mgq => 0,
         },
         rods => {
            uniqueness => { type => 'UQN', extension => '.uq'},
         },
      },
);

sub help_summary {
   return 'genotypes a set of samples';
}

sub help_text {
   my $self = shift;
   my $prog = $self->cl_name;
   my $summary = $self->help_summary;
   return <<EOM

Command 

  genotype - $summary

Syntaxis

  $prog genotye -g genotyper -i input1.bam -i input2.bam [-r REGION1 -r REGION2] -o output [further-genotyper-options]  


Genotypers

  Valid genotypers include:

  Samtools (mpileup/bcftools/vcfutils.pl)

  SamtoolsLegacy (pileup/varFilter)

  GATK 

EOM
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $config = Maltools::Config->new();
  my $output = "-"; my $errput = "-";
  my @regions = ();
  my $fragments = 1;
  my $genotyper = "DiploidGenotyper";
  my $bq_output = undef;
  my $baq_bq_output = undef;
  my $rbsq_output = undef;
  my $rmq_output = undef;
  my $reference = undef;
  my $max_depth = 100;
  my @files = ();
  my %options = ();
  my %rods = ();
  GetOptions(
    "input|i=s@" => \@files,
    "genotyper|G=s" => \$genotyper, 
    "t=i" => \$fragments, 
    "o=s" => \$output,
    "e=s" => \$errput, 
    "l=s@" => \@regions,
    "reference|ref|r=s" => \$reference,
    "rod=s" => \%rods,
    'option|opt=s%' => \%options);
  exists($genotypers{$genotyper}) or return $self->error_return("unsupported genotyper '$genotyper' ");
  
  $output ne "-" or return $self->error_return("you need to specify and output file name");
  $reference ||= $config->get_reference()->sequence_file;

  if ($genotyper eq 'Samtools') {
    return $self->samtools_genotyping(files => \@files, output => $output, reference => $reference, max_depth => $max_depth);
  }
  elsif ($genotyper eq 'SamtoolsLegacy') {
    return $self->samtools_legacy(files => \@files, output => $output, reference => $reference);
  }
  my $dc = data_config();
  $genotyper = $genotypers{$genotyper};
  my $tool = Maltools::Tool->tool_by_name($genotyper->{tool});
  my %genotyper_options = %{$genotyper->{options}};
  $genotyper_options{$_} = $options{$_} foreach  (keys %options);
  my $rods = $genotyper->{rods};
  if (keys %rods) {
     foreach my $k (keys %rods) {
        my $file = $rods{$k};
        my $type = "";
        $type = $1 if $k =~ s/,(\S+)$//;
        $rods->{$k} = { type => $type, file => $file };
     }
  }
  my $annotations = $genotyper->{annotations};
  if (exists($options{annotations}) || exists($options{annotation}) || exists($options{A})) {
    my @anno = ();
    push @anno, ( split(/\s*,\s*/,delete($options{annotation}) ) ) if  exists($options{annotation});
    push @anno, ( split(/\s*,\s*/,delete($options{annotations}) ) ) if  exists($options{annotations});
    push @anno, ( split(/\s*,\s*/,delete($options{A}) ) ) if  exists($options{A});
    $annotations = [ @anno ];
  }
  
  @ARGV = grep { $_ ne '--' } @ARGV;
  my $job = $tool->job(
      inputs => { 
          reference => $reference,
          genotyper => $genotyper->{genotyper},
          annotations => $annotations,
          rods => $genotyper->{rods},
          options => \%options,
          samples => [@files] , 
          @regions ? (intervals => [@regions]) : (),
          ($fragments > 1) ? ( threads => $fragments) : (),
          arguments => [@ARGV] }, 
      outputs => { 
          out => $output  });

  my $engine = $self->resolve_engine();
  print STDERR $self->engine_name,"\n";
  unless ($engine->run_job($job)) {
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotyping with message: " . $job->error_message);
  }
  else {
    return $self->ok_return;
  }
}

sub sample_to_file {
  my ($self,$file,$class,$dc) = @_;
  if ($class eq "file") {
    return undef unless -f $file;
    return undef unless -f $file . ".bai";
  }
  return $dc->sample_bam(ox_code => $file, indexed => 1, class => $class);
}

sub samtools_genotyping {
   my $self = shift;
   my %args = @_;
   my @bad = grep { ! -e $_ } @{$args{files}};
   return $self->error_return("some input files are not found: " . join(", ",@bad)) if $#bad >= 0;
   @bad = grep { ! -f $_ } @{$args{files}};
   return $self->error_return("some input files are not regular files: " . join(", ",@bad)) if $#bad >= 0;
   @bad = grep { ! -f $_ } @{$args{files}};
   return $self->error_return("some input files are not readable: " . join(", ",@bad)) if $#bad >= 0;

  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('samtools-genotyper');
  my $job = $tool->job(
        inputs => { in => $args{files}, max_depth => $args{max_depth} || 100, reference => $args{reference}, arguments => $args{others} },
        outputs => { out =>  $args{output} });
  unless ($engine->run_job($job)) {
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotping with message: " . $job->error_message);
  }
  else {
    return $self->ok_return;
  }
}

sub samtools_legacy {
   my $self = shift;
   my %args = @_;
   my @bad = grep { ! -e $_ } @{$args{files}};
   return $self->error_return("some input files are not found: " . join(", ",@bad)) if $#bad >= 0;
   @bad = grep { ! -f $_ } @{$args{files}};
   return $self->error_return("some input files are not regular files: " . join(", ",@bad)) if $#bad >= 0;
   @bad = grep { ! -f $_ } @{$args{files}};
   return $self->error_return("some input files are not readable: " . join(", ",@bad)) if $#bad >= 0;

  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('samtools-legacyGenotyper');
  my $job = $tool->job(
        inputs => { in => $args{files}, reference => $args{reference} },
        outputs => { out =>  $args{output} });
  unless ($engine->run_job($job)) {
    return $self->error_return("error occurred (code " . $job->return_code . ") during genotping with message: " . $job->error_message);
  }
  else {
    return $self->ok_return;
  }
}


1;
