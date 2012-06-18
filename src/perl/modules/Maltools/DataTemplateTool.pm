package Maltools::DataTemplateTool;

use base 'Maltools::Tool';
use strict;
use warnings;
use Text::Template;

our $OUTPUTS = { 
  out_dir => { type => 'directory', mandatory => 1 } # where to generat the classification information.
};

sub script_template {
   no strict 'refs';
   my $self = shift;
   my $class = ref($self);
   my $data = \*{$class . "::DATA"};
   return Text::Template->new (TYPE => 'FILEHANDLE', SOURCE => $data);
}


1;

