package Maltools::ToolData::flag;

use strict;
use warnings;

use Moose;

extends 'Maltools::ToolData::Type';

has '+name' => ( default => 'flag');

1;
