package Maltools::Reference::Location;

use Moose;

has reference => (is => 'ro', isa => 'Maltools::Reference', required => 1);
has sequence => (is => 'ro', isa => 'Str', required => 1);
has position => (is => 'ro', isa => 'Int', required => 1);

sub compare_with {
  my $self = shift;
  my $other = shift;
  $other->isa('Maltools::Reference::Location') 
    or die "the other need to be a location as well";
  my $reference = $self->reference;
#  $reference == $other->reference or die "cannot compare with a location from a different reference";
  my $my_chr = $self->sequence;
  my $ot_chr = $other->sequence;
  if ($other->sequence eq $self->sequence) {
    return $self->position <=> $other->position;
  }
  else {
    my $dictionary = $reference->dictionary;
    foreach my $entry (@$dictionary) {
       my $dict_chr = $$entry{name};
       return -1 if $dict_chr eq $my_chr;
       return 1 if $dict_chr eq $ot_chr;
    }
    die "did not found $my_chr nor $ot_chr in reference dictionary";
  }  
}

use overload 
   fallback => 1,
   "cmp" => sub { $_[2] ? $_[1]->compare_with($_[0]) : $_[0]->compare_with($_[1])  },
   '""'  => sub { sprintf('{%s %s}',$_[0]->sequence,$_[0]->position) };
1;
