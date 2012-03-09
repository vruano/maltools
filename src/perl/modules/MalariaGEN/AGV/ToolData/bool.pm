package MalariaGEN::AGV::ToolData::bool;

use strict;
use warnings;

use Moose;

extends 'MalariaGEN::AGV::ToolData::num';

has '+name' => ( default => 'bool');

1;
