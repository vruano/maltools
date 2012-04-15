package MalariaGEN::AGV::Commands::Concat;

use strict;
use warnings;
use base 'MalariaGEN::AGV::Command';

use Cwd qw(realpath);
use File::Basename qw(dirname);
use Getopt::Long;
use IO::File;

use MalariaGEN::AGV::Config;
use MalariaGEN::AGV::Tool;

sub help_summary {
   return 'concatenate files of the same type (currently only variant calls supported) into one';
}

sub help_text {
   my $program = $_[0]->cl_name;
   my $__or___ = '  or' . (' ' x (length($program) - 4));
   my $summary = $_[0]->help_summary;
  
   return <<EOM
Summary:
 
  concat - $summary

Synopsis:

  $program concat --vcf input1.vcf --vcf input2.vcf ... [-r reference.fa] -o output.vcf
  $__or___ concat vcf -i input1.vcf -i input2.vcf ... [-r reference.fa]  -o output.vcf
  $__or___ concat vcf -f input-files.list ... [-r reference.fa] -o output.vcf
  $__or___ concat --vcf-files inputs-files.list ... [-r reference.fa] -o output.vcf

Description:

  This command combine several files of some kind into a single output file. 

  For convencience you can indicate the input file type in two ways as indicated in 
  the Synopsis above. 

  If no reference is provided it will try to use the default reference based in
  the configuration. 

Options:

  <type> : (vcf) 

     indicates the type of the following input files that do not have a type 
     explicitly stated

  -i <file-name>
   
     input file without explicit type; this is determined by last <type> or
     file name extension otherwise.

  -o <file-name>

     output file name.

  -r <file-name>

     reference file name.

  --vcf <file-name>

     input vcf variation call file.

  --input-files <file-name>
    
     indicate a file that contains a list of inputs without specific type;
     this is determined by last <type> or each individual file name extension.

  --vcf-files <file-name>

     indicates a file that contains a list of input vcf variation call files.


EOM
;
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my @inputs = ();
  my $output = undef;
  my $type = '';
  my $reference = undef;
  my @input_lists = ();
  GetOptions(
     "<>" => sub { "$_[0]" },
     "input|i=s@" => sub { push @inputs, [ $type, $_[1] ] }, 
     "output|o=s" => \$output, 
     "vcf=s" => sub { push @inputs, [ 'vcf', $_[1] ] },
     "vcf-files=s@" => sub { push @input_lists, [ 'vcf', $_[1] ] },
     "input-files|f=s@" => sub { push @input_lists, [ $type, $_[1] ] }, 
     "ref|r=s" => \$reference, );
  
  $self->_process_input_lists(\@input_lists,\@inputs);
  $@ and return $self->error_return($@);
  
  @inputs = map {
    my ($type,$file) = @$_;
    unless ($type) {
    	$type = $1 if $file =~ /\.([^\.]+)$/;
    }
    [ lc($type), $file ];
  } @inputs;


  return $self->error_return("no input was specified") unless @inputs;

  return $self->error_return("no output was specified") unless $output;

  return $self->error_return("output exists and is not a regular file") if -e $output && ! -f $output;

  return $self->error_return("output's parent directory does not exists or is unreachable") unless -e dirname($output);

  return $self->error_return("output's parent directory is in fact not a directory") unless -d dirname($output);

  return $self->error_return("reference file does not exists") if $reference && !-e $reference;
 
  return $self->error_return("reference file is not a regular file") if $reference && -e $reference && !-f $reference;
 
  $reference = $reference ? MalariaGEN::AGV::Reference->new(sequence_file => $reference) : MalariaGEN::AGV::Config->new()->get_reference();
  
  $reference or return $self->error_return("there is no default reference, you must indicate one");
  $reference->exists() or return $self->error_return("reference sequence file could not be found '$reference'");

  my @bad_types = grep { $_ && $_ ne 'vcf' } map { $$_[0] } @inputs;
  return $self->error_return("unknown input type '$bad_types[0]'") if @bad_types;
 
  my @unknown_type_files = map { $$_[1] } grep { ! $$_[0] } @inputs;
  return $self->error_return("cannot determine type for input files '" . join(", ", @unknown_type_files)."'") if @unknown_type_files;

  my @non_existent_files = grep { ! -e $_ } map { $$_[1] } @inputs;
  return $self->error_return("some input files cannot be reached or do not exists '" . join(", ", @non_existent_files) . "'") if @non_existent_files;

  my @non_file_inputs = grep { ! -f $_ } map { $$_[1] } @inputs;
  return $self->error_return("some inputs are not regular files '" . join(", ", @non_file_inputs) . "'") if @non_file_inputs;

  my %inputs_per_type = ();

  push @{$inputs_per_type{$$_[0]}}, $$_[1] foreach @inputs;

  return $self->error_return("mixing incompatible input types '" . join(", ", keys %inputs_per_type) . "'") if (keys %inputs_per_type > 1); 
  
  $type = (keys %inputs_per_type)[0];

  my $job;
  if ($type eq 'vcf') {
    $job = $self->concat_vcf_job(inputs => $inputs_per_type{vcf}, output => $output, reference => $reference);
  }
  else { # must never happen as is detected before but just in case....
    return $self->error_return("unknown input type '$type'");
  }
  push @inputs,@ARGV;
  
  my $engine = $self->resolve_engine;
  $engine->run_job($job) 
	or $self->error_return("error executing the align command with exit code : " . $job->return_code,$job->error_message );
  return $self->ok_return;
}

sub _process_input_lists {
   my ($self,$lists,$inputs) = @_;

   foreach my $l (@$lists) {
      my ($type,$file_name) = @$l;
      my $num = 0;
      print STDERR "processing '$file_name' ...\n";
      open my $fh , $file_name or return $@ = "could not open '$file_name' for reading";
      while (my $i_file_name = <$fh>) {
        chomp $i_file_name;
        next unless $i_file_name =~ /\S/;
        next if $i_file_name =~ /^#/;   
        $num++;
        push @$inputs , [ $type, $i_file_name ];
      }
      close $fh; 
      print STDERR "warning! '$file_name' does not seem to contain any input file name!\n" unless $num;
   }
}

sub concat_vcf_job {
  my $self = shift;
  my %opts = @_;
  my $tool = MalariaGEN::AGV::Tool->tool_by_name('vcf-concat');
  my $job = $tool->job(inputs => { reference => $opts{reference}, inputs => $opts{inputs} }, outputs => { out => $opts{output}});
  return $job;
}








1;
