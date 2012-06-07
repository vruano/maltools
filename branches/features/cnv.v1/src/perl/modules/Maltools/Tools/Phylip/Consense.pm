package Maltools::Tools::Phylip::Consense;

use strict;
use warnings;
use Moose;

extends 'Maltools::Tool';

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift, @_) };

our $INPUTS = {
   in => { type => 'file', mandatory => 1},
};

our $OUTPUTS = {
   out_tree => { type => 'file', default => "-" },
   out_file => { type => 'file', default => "-" },
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

my $input = "{ $J->input("in") }";
my $outfile = "{$J->output("out_file")}";
my $outtree = "{$J->output("out_tree")}";

my $wd = tempdir("consense_XXXX",DIR => getcwd);
my $intree = catfile($wd,"intree");
my $output_outfile = catfile($wd,"outfile");
my $output_outtree = catfile($wd,"outtree");

sub resolve_intree \{
   if ($input eq "-") \{
     my $fh = IO::File->new($intree,"w") or die "could not open temporal intree file '$intree'";
     while (my $line = <STDIN>) \{
        print $fh $line;
     \}
     $fh->close();
   \}
   else \{
     return if link($input,$intree); # A hard link worked well. 
     return if symlink($input,$intree); # A symbolic link worked well.
     copy($input,$intree) or die "could not copy input '$input' as the intree file '$intree'";
   \}
\}


## Main code;

resolve_intree();
print STDERR "Processing input trees ...";


my $command = "cd $wd; consense 2>&1 > /dev/null <<EOF\n";
  $command .= "2\n"; # remove progress indication in output
  $command .= "Y\n";
  $command .= "EOF\n";

system($command);
if ($?) \{
  print STDERR "error\n";
  my $mess = "problems executing the consense command with return code '$?' and error message '$@'";
  #eval \{ remove_tree($wd) \};
  die $mess;
\}
else \{
  print STDERR "done\n";
  # copy outputs.
  if ($outfile) \{
    if ($outfile eq "-") \{
      copy($output_outfile,\*STDOUT) or die "could not dump outfile '$output_outfile' into the standard output";
    \}
    else \{
      move($output_outfile,$outfile) or die "could not move outfile '$output_outfile' into requested location '$outfile'";
    \}
  \}
  if ($outtree) \{
    if ($outtree eq "-") \{
      copy($output_outtree,\*STDOUT) or die "could not dump outtree '$output_outtree' into the standard output";
    \}
    else \{
      move($output_outtree,$outtree) or die "could not move outtree '$output_outtree' into requested location '$outtree'";
    \}
  \}
\}
# won't remove files if there is any error.
remove_tree($wd);
