package MalariaGEN::AGV::Tools::Copy;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift) }
  
my $INPUTS = {
    src => {  type => 'file', mandatory => 1 },
};

my $OUTPUTS = { 
    dest => { type => 'file', mandatory => 1}
};

sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
   return 'sh'; 
}

sub script_template {
   return TTS('cp { $J->input("src") } { $J->output("dest") }' );
}

1;
