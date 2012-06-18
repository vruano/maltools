package Maltools::Reference::ROD;

use Moose;
use File::Spec::Functions qw(catfile);

has 'file' => (is => 'ro', isa => 'Str', required => 1);
has 'type' => (is => 'ro', isa => 'Str', predicate => 'has_type');

1;
