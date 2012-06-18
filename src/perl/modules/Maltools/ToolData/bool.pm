package Maltools::ToolData::bool;

use strict;
use warnings;

use Moose;

extends 'Maltools::ToolData::num';

has '+name' => ( default => 'bool');

1;
