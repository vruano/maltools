package Samtrak::Data::File::Manager;

use strict;
use warnings;
use Samtrak::Data::File;
use Samtrak::DB::Manager;

use base 'Samtrak::DB::Manager';

sub object_class {'Samtrak::Data::File'};

__PACKAGE__->make_manager_methods('files');

1;
