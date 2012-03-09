package MalariaGEN::AGV::Tools::Bwa::Aln;

use base 'MalariaGEN::AGV::Tool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING', shift, @_) }
  
my $INPUTS = {
    query => {  type => 'file', mandatory => 1 },
    db => {  type => 'file', mandatory => 1 },
};

my $OUTPUTS = { 
    out => { type => 'file', mandatory => 1,
             default => TTS('{my $q = $J->input("query"); 
                                 $q =~ s/(\.fa)|(\.fasta)$//; 
                                 return $q .= ".sai"}') },
};

sub new {
  my ($class,%args) = @_;
  $args{inputs} = $INPUTS unless defined $args{inputs};
  $args{outputs} = $OUTPUTS unless defined $args{outputs};
  return $class->SUPER::new(%args);
}

sub interpreter {
   my ($self,%args) = @_;
   return TTS('{$SH}'); 
}

sub script {
   my ($self,%args) = @_;
   return TTS('bwa aln { $J->input("db") } { $J->input("query") } > { $J->output("out") }' );
}

1;
