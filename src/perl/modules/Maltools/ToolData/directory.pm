package Maltools::ToolData::directory;

use strict;
use warnings;
use Moose;

extends 'Maltools::ToolData::file';

has '+has_files' => ( default => 1 );

1;
