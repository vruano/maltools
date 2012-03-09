package MalariaGEN::AGV::Commands::SampleTable;

use Moose;
extends 'MalariaGEN::AGV::Command';

use strict;
use warnings;

use MalariaGEN::AGV::Config qw(sanger_config);
use Getopt::Long qw(:config pass_through);

use IO::File;
use File::Spec::Functions qw(catfile);
use base 'MalariaGEN::AGV::Command';


has 'data_config' => (is => 'ro', lazy => 1 , default => sub { MalariaGEN::AGV::Config::data_config() });

sub help_summary {
   return 'generates a sample table from a VCF giving general stats per sample including its classificiation between GOOD and POOR';
}

sub help_text {
   return "sample-table:\n\t" .
          "execute some query relevant to the AGV pipeline\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " sample-table [-i input.vcf | < input.vcf ] [-o output.tsv | > output.tsv ]\n";
}

sub execute {
    my ($self,$samtrak,$params) = @_;
    local @ARGV;
    @ARGV = @{$params->{arguments}};
    my $output = "-";
    my $input = "-";
    GetOptions("o=s" => \$output, "i=s" => \$input);
    return $self->error_return("could not find the input file '$input'") if ($input ne "-" && !-e $input);
    my $engine = $self->resolve_engine;
    my $tool = MalariaGEN::AGV::Tool->tool_by_name('vcf-goodSamples') or return $self->error_return("good sample calculation error","could not find a tool with id 'vcf-goodSamples'");
    my $job = $tool->job(inputs => { in => $input }, outputs => { out => $output } );
    unless ($engine->run_job($job)) {
       return $self->error_return("error executing vcf-goodSamples tool with return code " . $job->return_code . " and message: " . ($job->error_message || 'no-message'));
    }
    return $self->ok_return;
}

1;
