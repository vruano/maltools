package MalariaGEN::AGV::Tools::Phylip::Neighbor;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
   in => { type => 'file', mandatory => 1},
   format => { type => 'string', default => 'lower-triangular' }, # lower-triangular, upper-triangular or squared.
   ds_count => { type => 'num' , default => undef }, # indicate how many datasets are found in 'in', leave it undef for the program to find-out itself.
   seed => { type => 'num', default => 13 },
};

our $OUTPUTS = {
   out_tree => { type => 'file', default => "-" },
   out_file => { type => 'file', default => "" },
};

sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
  return "perl";
}


sub command_template {
   return TTS('{$T->interpreter} {$S}');
}


sub script_template {
  return Text::Template->new(TYPE=> 'FILEHANDLE', SOURCE => \*DATA);
}


1;

__DATA__

use strict;
use warnings;

use Cwd qw(getcwd);
use File::Temp qw(tempdir);
use IO::File;
use File::Copy qw(copy move);
use File::Path qw(remove_tree);
use File::Spec::Functions qw(catfile);

my $input = "{$J->input("in")}";
my $input_format = "{$J->input("format")}";
my $dataset_count = undef; # undef means that you need to guess it.
my $outfile = "{$J->output("out_file")}";
my $outtree = "{$J->output("out_tree")}";
my $seed = {$J->input("seed") || 13};

my $wd = tempdir("neighbor_XXXX",DIR => getcwd);
my $infile = catfile($wd,"infile");
my $neighbor_of = catfile($wd,"outfile");
my $neighbor_ot = catfile($wd,"outtree");

sub resolve_infile \{
   if ($input eq "-") \{
     my $fh = IO::File->new($infile,"w") or die "could not open temporal infile '$infile'";
     my $ds_count = 0;
     while (my $line = <STDIN>) \{
        $ds_count++ if ($line =~ /^\s*(\d+)\s*$/);
        print $fh $line;
     \}
     $fh->close();
     if (defined $dataset_count && $dataset_count != $ds_count) \{
       die "there is a discrepancy between the declared number of datasets '$dataset_count' and the actual number '$ds_count'";
     \}
     else \{
       $dataset_count = $ds_count;
     \}
   \}
   else \{
     return if link($input,$infile); # A hard link worked well. 
     return if symlink($input,$infile); # A symbolic link worked well.
     copy($input,$infile) or die "could not copy input '$input' as the infile '$infile'";
   \}
\}



sub resolve_dataset_count \{
   return $dataset_count if ($dataset_count);
   # input must be a file and so infile must have the same value.
   open(INFILE,$infile) or die  "could not open input file '$infile'";
   my $ds_count = 0;
   while (my $line = <INFILE>) \{
     $ds_count++ if ($line =~ /^\s*(\d+)\s*$/);
   \}
   close(INFILE);
   return $dataset_count = $ds_count;
\}

## Main code;

resolve_infile();
resolve_dataset_count();
print STDERR "Processing $dataset_count datasets ...";

my $command = "cd $wd; neighbor 2>&1 > /dev/null <<EOF\n";
  $command .= "L\n" if $input_format eq "lower-triangular";
  $command .= "R\n" if $input_format eq "upper-triangular";
  $command .= "M\n";
  $command .= "$dataset_count\n";
  $command .= "$seed\n";
  $command .= "2\n"; # remove progress indication in output
  $command .= "Y\n";
  $command .= "EOF\n";

system($command);
if ($?) \{
  print STDERR "error\n";
  my $mess = "problems executing the neighbor command with return code '$?' and error message '$@'";
  #eval \{ remove_tree($wd) \};
  die $mess;
\}
else \{
  print STDERR "done\n";
  # copy outputs.
  if ($outfile) \{
    if ($outfile eq "-") \{
      copy($neighbor_of,\*STDOUT) or die "could not dump outfile '$neighbor_of' into the standard output";
    \}
    else \{
      move($neighbor_of,$outfile) or die "could not move outfile '$neighbor_of' into requested location '$outfile'";
    \}
  \}
  if ($outtree) \{
    if ($outtree eq "-") \{
      copy($neighbor_ot,\*STDOUT) or die "could not dump outtree '$neighbor_ot' into the standard output";
    \}
    else \{
      move($neighbor_ot,$outtree) or die "could not move outtree '$neighbor_ot' into requested location '$outtree'";
    \}
  \}
\}
# won't remove files if there is any error.
remove_tree($wd);
