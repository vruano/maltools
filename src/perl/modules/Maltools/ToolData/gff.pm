package Maltools::ToolData::gff;

use strict;
use warnings;
use Moose;

extends 'Maltools::ToolData::rod';

has '+has_files' => ( default => 1 );

1;
