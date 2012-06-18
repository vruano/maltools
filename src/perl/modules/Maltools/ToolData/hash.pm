package Maltools::ToolData::hash;

use strict;
use warnings;

use Moose;

extends 'Maltools::ToolData::Type';

has '+name' => ( default => 'hash');

1;
