package Maltools::Tools::CrossCoverage;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use vars qw(%ENV);
use File::Basename qw(dirname);
use Text::Template;
use IO::File;
use Cwd qw(realpath);
use JSON::XS;
use File::Spec::Functions qw(catfile file_name_is_absolute);
use POSIX;

our $INPUTS = {
   in_files => { type => 'file', multiple=> 1, mandatory => 1 }
};

our $OUTPUTS = { 
   out_dir => { type => 'file', mandatory => 1 }
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 8000;
}

sub calculate_cpu_time {
  return 60 * 60 * 1;
}

sub interpreter {
  return 'Rscript';
}

1;

__DATA__
{ 
  @in_files = @{$J->input('in_files')};
  $out_dir = $J->output('out_dir'); '' }

list.files = c("{join("\",\"",@in_files)}");

num.files = length(list.files)

for (i in 1:num.files)
\{
	path.to.file = list.files[i];
	cov = read.table(path.to.file,skip = 1, header = TRUE, sep ="\t");

	if (i == 1)
	\{	row.sums = rowSums(cov[,4:7]);	\}
	else
	\{	row.sums = row.sums + rowSums(cov[,4:7]); \}
\}

row.sums = row.sums/num.files;

cov.quants =quantile(row.sums,probs = seq(0.01,1,by=0.02));
write.table(file = "{$out_dir}/cov.quants", cov.quants, row.names = TRUE, col.names = FALSE, quote = FALSE, sep ="\t");

for (i in 1:14)
\{
	chromo = paste("MAL",i,sep="");
	this.chromo = which(cov[,1]==chromo);
	file.out = paste("{$out_dir}",paste("coverage.MAL",i,sep=""),sep="/");

	write.table(cbind(cov[this.chromo,1:2],row.sums[this.chromo]),file = file.out, 
          sep="\t",col.names = c("Chromosome","Position", "Coverage"), 
          row.names= FALSE, quote= FALSE);
\}


