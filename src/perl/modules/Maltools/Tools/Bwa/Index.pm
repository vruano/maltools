package Maltools::Tools::Bwa::Index;

use base 'Maltools::Tool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING',shift,@_); }

sub new {
  my ($class,%args) = @_;
  $args{inputs} = {
    alg => {  type => 'enum' , mandatory => 0, values => [qw[is bwtsw]], default => 'bwtsw' },
    cs_idx => {  type => 'bool', mandatory => 0, default => 0 },
    in => {  type => 'file', mandatory => 1 },
  } unless defined $args{inputs};
  $args{outputs} = {
    out => { type => 'fileset', mandatory => 0, 
             default => Text::Template->new(TYPE => 'STRING', SOURCE => '{$J->input("in")}') },
  } unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
  return TTS('{$SH}');
}

sub script {
   my ($self,%args) = @_;
   return TTS('bwa index { $J->has_output("out") ? "-p " . $J->output("out") : "" }' .
                  '-a { $J->input("alg") } { $J->input("cs_idx") ? "-c" : ""} ' . 
                  '{ $J->input("in") }');  
}

1;
