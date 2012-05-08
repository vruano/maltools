package Maltools::ToolData::string;

use strict;
use warnings;

use Moose;

extends 'Maltools::ToolData::Type';

has '+name' => ( default => 'string');

1;
