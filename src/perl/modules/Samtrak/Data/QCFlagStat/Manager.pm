package Samtrak::Data::QCFlagStat::Manager;

use strict;
use warnings;

use base 'Samtrak::DB::Manager';


sub object_class {'Samtrak::Data::QCFlagStat'};

__PACKAGE__->make_manager_methods('qcflagstats');

sub get_by_oxford_code {
   my ($self,$ox_code,@options) = @_;
   my $results = $self->get_sample2lanes(query => [ ox_code => $ox_code ],@options);
   return [ ] unless defined $results;
   return [ $results ] if ref($results) =~ /Samtrak/;
   return $results;
}

1;
