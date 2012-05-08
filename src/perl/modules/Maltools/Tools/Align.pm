package Maltools::Tools::Align;

use base 'Maltools::DataTemplateTool';
use strict;
use warnings;
use Text::Template;

sub TTS { Text::Template->new(TYPE => 'STRING', SOURCE => shift) }
  
my $INPUTS = {
    in => {  type => 'bam' , multiple => 1, mandatory => 1 },
    force => { type => 'bool', default => 0},
    ref => { type => 'reference', mandatory => 1}
};

my $OUTPUTS = { 
    out => { type => { name => 'bam', indexed => 1} , mandatory => 1}
};

sub interpreter {
   return 'perl'; 
}


1;
__DATA__
use strict;
use warnings;

use Maltools::Alignment::BWA qw(bwa);

bwa(in => "{$J->input('in')}",
    out => "{$J->output('out')}",
    ref => "{$J->input('ref')}",
    sort => 1,
    index => 1);
