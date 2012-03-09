package MalariaGEN::AGV::Tools::Bwa::Sampe;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING',shift,@_); }

my $INPUTS = {
  db => { type => 'file', mandatory => 1 },
  f_sai => { type => 'file', mandatory => 1 },
  r_sai => { type => 'file', mandatory => 1 },
  f_fq => { type => 'file', mandatory => 1},
  r_fq => { type => 'file', mandatory => 1},
};

my $OUTPUTS = {
  out => { type => 'file', mandatory => 1 }
};

sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
  return TTS('{$SH}');
}

sub script {
   my ($self,%args) = @_;
   return TTS('bwa sampe { $J->input("db") } { $J->input("f_sai") } { $J->input("r_sai") } { $J->input("f_fq") } { $J->input("r_fq") } > { $J->output("out") }');
}

1;
