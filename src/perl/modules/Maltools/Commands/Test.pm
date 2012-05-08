package Maltools::Commands::Test;

use Moose;

extends 'Maltools::Command';

use strict;
use warnings;

use Maltools::Config qw(data_config 
       sanger_config reference_config jobs_config);
use Maltools::Tool;
use Getopt::Long;
use IO::File;

has '+engine_name' => ( default => 'local' );


sub hidden {
   return 1;
}

sub help_summary {
   return 'test command gives out the engine name';
}

sub help_text {
   return "test\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " test\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $copy = Maltools::Tool->tool_by_name("test");
  my $jc = jobs_config();
  my $engine = $jc->engine();
  print "Default engine: " . $engine->name,"\n";
  return $self->ok_return;
}


1;
