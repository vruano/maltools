package MalariaGEN::AGV::ToolData::gff;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::ToolData::rod';

has '+has_files' => ( default => 1 );

1;
