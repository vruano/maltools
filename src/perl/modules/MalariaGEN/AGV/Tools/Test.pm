package MalariaGEN::AGV::Tools::Test;

use strict;
use warnings;
use Text::Template;

use Moose;

extends 'MalariaGEN::AGV::Tool';

our $OUTPUTS = {
  out => { type => 'file', default => "-" },
};

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift) }
  
sub interpreter {
   return 'sh'; 
}

sub script_template {
   return TTS('{$out = $J->output("out") || "-"; $out = "/dev/stdout" if $out eq "-"; "" }echo The engine is { $E->name } > {$out}' );
}

1;
