package MooseX::Cloneable;

use Moose::Role;

sub clone {
    my $self = shift;
    my %params = @_;
    return $self->meta->clone($self,%params);
}

1;