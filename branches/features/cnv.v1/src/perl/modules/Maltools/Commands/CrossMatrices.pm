package Maltools::Commands::CrossMatrices;

use strict;
use warnings;

use Moose;
use File::Copy qw(copy);
use Cwd qw(realpath);

use Maltools::Config qw(data_config 
  sanger_config reference_config jobs_config);

use Getopt::Long;
extends 'Maltools::Command';


sub hidden {
  return 1;
}

sub help_summary {
   return 'constructs cross matrices directory from a given variants file (VCF)';
}

sub help_text {
   return "snps-filter\n\t" .
      $_[0]->help_summary . "\n" .
      "Syntaxis:\n" .
      $_[0]->cl_name . " cross-matrices -i input.vcf -o out-dir [-p PARENT]{2,} [--include all|FILTER]* [--exclude all|FILTER]*\n" .
      "\t\t--include indicates ''filters'' that should be included in the output. 'all' means that all filters should be included by default\n" .
      "\t\t--exclude indicates ''filters'' that should be excluded from the output. 'all' means that all filters should be excluded by default\n" .
      "\t\tNotice that PASS LIT CREDIBLE AND TYPABLE are also valid values for include or exclude\n" .
      "\t\t--excl-invs indicates whether sites that do not present variability across samples should be excluded from the output\n".
      "\t\t\tno : no exclussion (default)\n".
      "\t\t\talleles : exclusion of all sites without variance at the allele call level\n".
      "\t\t\tgenotypes : exclusion of all sites without variance at the genotype calls\n". 
      "Examples:\n" .
      "\t" . $_[0]->cl_name . " cross-matrices -i input.vcf -o out-dir [--exclude all]? --include PASS\n" .
      "\t\t\tIncludes only 'PASS' variants. --exclude all is implicit when using any 'include'\n" .
      "\t" . $_[0]->cl_name . " cross-matrices -i input.vcf -o out-dir [--exclude all]? --include CREDIBLE\n" .
      "\t\t\tIncludes only credible variants (with a 'PASS', with credible or typable or lit info)\n" .
      "\t" . $_[0]->cl_name . " cross-matrices -i input.vcf -o out-dir --include MinCoverage\n" .
      "\t\t\tIncludes only variants that fail the min. coverage filter\n" ;
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $manifest = undef;
  my $rc = reference_config();
  my $input = undef;
  my $output = undef;
  my $excl_invs = 'no';
  my @includes = ();
  my @excludes = ();
  my @parents = ();
  GetOptions("input|i=s" => \$input, "excl-invs|ei=s" => \$excl_invs, "output|o=s" => \$output, 
             "parent|p=s" => \@parents, "include|inc=s@" => \@includes, "exclude|excl=s" => \@excludes);

  $excl_invs eq "no" || $excl_invs eq "alleles" || $excl_invs eq "genotypes" or return $self->error_return("illegal value for --excl-invs '$excl_invs'");
  -f $input or return $self->error_return("you must provide an input variant file");
  -e $output && ! -d $output and return $self->error_return("the ouput directiory exists already but is not a directory");
  scalar(@parents) >= 2 or return $self->error_return("there must be at least two parents"); 
   
   if ($#includes < 0 && $#excludes < 0) {
     push @includes, 'CREDIBLE';
     push @excludes, 'ALL';
   } 
   elsif ($#includes >= 0 && $#excludes < 0) {
     push @excludes, 'ALL' unless grep { uc($_) eq "ALL" } @includes;
   }
   elsif ($#includes < 0 && $#excludes >= 0) {
     push @includes, 'ALL' unless grep { uc($_) eq "ALL" } @excludes;
   }

   my %includes = map { uc($_) => 1 } @includes;
   my %excludes = map { uc($_) => 1 } @excludes;

   my @collitions = grep { exists($excludes{$_}) } keys %includes;

   @collitions and return $self->error_return("inclussions and exclussions cannot make reference to the same variant qualifiers: " . join(", ",@collitions));
  
  my $engine = $self->resolve_engine();
  my $tool = Maltools::Tool->tool_by_name('varToCrosses');
  my $job = $tool->job(inputs => { in => $input, invariants => $excl_invs, 
         parents => \@parents, includes => \@includes, excludes => \@excludes} , outputs => { out_dir => $output });
  my $result = $engine->run_job($job);
  unless ($result) {
    return $self->error_return("error occurred during cross matrices building",$job->error_message());
  }
  return $self->ok_return;
}


1;
