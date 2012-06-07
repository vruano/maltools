package Maltools::ToolData::num;

use strict;
use warnings;

use Moose;

extends 'Maltools::ToolData::string';

has '+name' => ( default => 'num');

1;
